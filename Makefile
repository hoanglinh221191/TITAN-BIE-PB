# ==========================================
# TITAN PROJECT MAKEFILE (Cross-platform)
# Windows: Clang (C++23) | Linux: GCC/Clang (C++23)
# Tích hợp Auto Dependency Tracking (.d)
# ==========================================

PROFILE ?= native
TOOLCHAIN ?= auto
TITAN_CXX23_FLOOR ?= clang20-gxx13
LLVM_VERSION ?= 22
LLVM_PREFIX ?= /usr/lib/llvm-$(LLVM_VERSION)
BUILD_ROOT = build
BUILD_DIR = $(BUILD_ROOT)/$(PROFILE)
SRC_ROOT = src
INCLUDE_DIRS = -I$(SRC_ROOT) -I$(SRC_ROOT)/gmxtraj/include -I$(SRC_ROOT)/eigen

COMMON_CXX_DEFINES = -DTITAN_CXX23_FLOOR_CLANG20_GXX13=1 -DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1
ifeq ($(PROFILE),native)
    COMMON_CXX_DEFINES += -DTITAN_BUILD_PROFILE_TOKEN=native
else ifeq ($(PROFILE),portable)
    COMMON_CXX_DEFINES += -DTITAN_BUILD_PROFILE_TOKEN=portable
else
    COMMON_CXX_DEFINES += -DTITAN_BUILD_PROFILE_TOKEN=unknown
endif

LOWMEM_CXXFLAGS = $(filter-out -O3 -pipe,$(CXXFLAGS)) -O1

ifeq ($(PROFILE),native)
    ARCH_FLAGS = -march=native
else ifeq ($(PROFILE),portable)
    ARCH_FLAGS = -mavx2 -mfma -mtune=generic
else
    $(error Unsupported PROFILE '$(PROFILE)'. Use PROFILE=native or PROFILE=portable)
endif

ifeq ($(OS),Windows_NT)
    # ─── CẤU HÌNH CHO WINDOWS (CLANG) ───
    MSYS2_UCRT64_BIN ?= C:/msys64/ucrt64/bin
    CXX  = $(MSYS2_UCRT64_BIN)/clang++.exe
    RM   = del /Q /F
    TARGET = titan.exe
    PIPE_FLAG = -pipe
    STACK_NOTE_FLAG =
    
    # Clang Flags: Thêm cờ chặn cảnh báo rác của Clang
    # Force Eigen to use its own aligned allocator on Windows.
    # This avoids AVX-aligned loads/stores from touching only 16-byte-aligned buffers.
    CXXFLAGS = $(INCLUDE_DIRS) -O3 $(ARCH_FLAGS) -fopenmp $(PIPE_FLAG) $(STACK_NOTE_FLAG) -std=c++23 -Wall -Wno-unused-command-line-argument -MMD -MP $(COMMON_CXX_DEFINES)
    
    # Ép dùng LLD linker siêu tốc của LLVM
    LDFLAGS = -fuse-ld=lld -static -static-libgcc -static-libstdc++ -fopenmp
    
    # Clang dùng ThinLTO đa luồng
    LTO_FLAG = -flto=thin
    
    # Lệnh dọn dẹp trên Windows CMD/MSYS2
    CLEAN_CMD = if exist "$(subst /,\,$(BUILD_DIR))" rmdir /S /Q "$(subst /,\,$(BUILD_DIR))" & if exist "$(TARGET)" del /Q /F "$(TARGET)"
    CLEAN_ALL_CMD = if exist "$(subst /,\,$(BUILD_ROOT))" rmdir /S /Q "$(subst /,\,$(BUILD_ROOT))" & if exist "$(TARGET)" del /Q /F "$(TARGET)"
    CLEAN_GENERATED_CMD = $(RM) network_edges.csv network_nodes.csv network_communities.csv network_community_summary.csv network_supergraph.csv pca_fel.csv tica_projections.csv extract_basins.sh tmp_test_dham*.in 2>nul || exit 0
else
    # ─── CẤU HÌNH CHO LINUX (GCC 13) ───
    TOOLCHAIN_NAME := $(TOOLCHAIN)
    ifeq ($(TOOLCHAIN_NAME),auto)
        ifneq ($(findstring clang,$(notdir $(CXX))),)
            TOOLCHAIN_NAME := clang
        else
            TOOLCHAIN_NAME := gcc
        endif
    endif
    ifeq ($(TOOLCHAIN_NAME),clang)
        ifeq ($(origin CXX),default)
            CXX = clang++-22
        endif
        ifeq ($(origin CC),default)
        endif
        OPENMP_FLAG = -fopenmp=libomp
        CXX_STDLIB_FLAG = -stdlib=libc++
        CLANG_STDLIB_INCLUDE_FLAGS = -nostdinc++ -isystem $(LLVM_PREFIX)/include/c++/v1 -isystem $(LLVM_PREFIX)/include/x86_64-pc-linux-gnu/c++/v1
        CLANG_OMP_INCLUDE_FLAGS = -isystem $(LLVM_PREFIX)/lib/clang/$(LLVM_VERSION)/include
        CLANG_STDLIB_LINK_FLAGS = -L$(LLVM_PREFIX)/lib -Wl,-rpath,$(LLVM_PREFIX)/lib
        CLANG_STDLIB_LIBS = -lc++abi
        STDLIB_LINK_FLAGS =
        LTO_FLAG = -flto=thin
    else ifeq ($(TOOLCHAIN_NAME),gcc)
        ifeq ($(origin CXX),default)
            CXX = g++
        endif
        ifeq ($(origin CC),default)
        endif
        OPENMP_FLAG = -fopenmp
        CXX_STDLIB_FLAG =
        CLANG_STDLIB_INCLUDE_FLAGS =
        CLANG_OMP_INCLUDE_FLAGS =
        CLANG_STDLIB_LINK_FLAGS =
        CLANG_STDLIB_LIBS =
        STDLIB_LINK_FLAGS = -static-libgcc -static-libstdc++
        LTO_FLAG = -flto=auto
    else
        $(error Unsupported TOOLCHAIN '$(TOOLCHAIN)'. Use TOOLCHAIN=auto, gcc, or clang)
    endif
    RM   = rm -f
    TARGET = titan
    PIPE_FLAG =
    STACK_NOTE_FLAG = -Wa,--noexecstack
    
    # GCC Flags: Chuẩn C++23, loại bỏ các cờ đặc thù của Clang
    # Keep Eigen on its own aligned-allocation path across Linux builds too.
    CXXFLAGS = $(INCLUDE_DIRS) -O3 $(ARCH_FLAGS) $(OPENMP_FLAG) $(CXX_STDLIB_FLAG) $(CLANG_STDLIB_INCLUDE_FLAGS) $(CLANG_OMP_INCLUDE_FLAGS) $(PIPE_FLAG) $(STACK_NOTE_FLAG) -std=c++23 -Wall -MMD -MP $(COMMON_CXX_DEFINES)
    
    # Linker mặc định của GCC (ld)
    LDFLAGS = $(STDLIB_LINK_FLAGS) $(CLANG_STDLIB_LINK_FLAGS) $(OPENMP_FLAG) $(CXX_STDLIB_FLAG) $(CLANG_STDLIB_LIBS)
    
    # GCC dùng LTO tự động nhận diện số nhân CPU
    
    # Lệnh dọn dẹp trên Linux bash
    CLEAN_CMD = rm -rf $(BUILD_DIR) $(TARGET)
    CLEAN_ALL_CMD = rm -rf $(BUILD_ROOT) $(TARGET)
    CLEAN_GENERATED_CMD = $(RM) network_edges.csv network_nodes.csv network_communities.csv network_community_summary.csv network_supergraph.csv pca_fel.csv tica_projections.csv extract_basins.sh tmp_test_dham*.in 2>/dev/null || true
endif

# ─── Tệp nguồn (Sources) ───
ROOT_SRCS = $(wildcard $(SRC_ROOT)/*.cpp)
PBE_GRID_SRCS = $(wildcard $(SRC_ROOT)/pbe/grid/*.cpp)
PBE_SRCS  = $(wildcard $(SRC_ROOT)/tabi/*.cpp) $(PBE_GRID_SRCS)
PBCTOPO_SRCS = $(wildcard $(SRC_ROOT)/pbctopo/*.cpp)
GMXTRAJ_SRCS = $(wildcard $(SRC_ROOT)/gmxtraj/*.cpp)
SRCS  = $(ROOT_SRCS) $(PBE_SRCS) $(PBCTOPO_SRCS) $(GMXTRAJ_SRCS)
OBJS  = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS  = $(OBJS:.o=.d)

# ==========================================
# BUILD RULES
# ==========================================

.DELETE_ON_ERROR:
.PHONY: all clean clean-all clean-generated release native portable release-native release-portable print-config

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo Linking $(TARGET) [$(PROFILE)]...
	@$(CXX) -o $@ $^ $(LDFLAGS)
	@echo Build successful!

$(BUILD_DIR)/$(SRC_ROOT)/titan_rf_model.o: CXXFLAGS := $(filter-out -O3,$(CXXFLAGS)) -O0

ifneq ($(OS),Windows_NT)
$(BUILD_DIR)/$(SRC_ROOT)/solv_pbc_pbctopo.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/anchor_graph.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/component_graph.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/component_refine.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/component_search.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/component_search_pair_table.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/llps_local.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pbctopo/local_graph.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/engine_postproc.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/impl_engine.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/solv_pbc_new_al.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/find_pull_direction.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/pull_direction_pathway.o: CXXFLAGS := $(LOWMEM_CXXFLAGS)
$(BUILD_DIR)/$(SRC_ROOT)/titan_rf_model.o: CXXFLAGS := $(filter-out -O3 -pipe,$(CXXFLAGS)) -O0
endif

$(BUILD_DIR)/%.o: %.cpp
	@echo Compiling $<...
ifeq ($(OS),Windows_NT)
	@if not exist "$(subst /,\,$(@D))" mkdir "$(subst /,\,$(@D))"
else
	@mkdir -p $(@D)
endif
	@$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo Cleaning up...
	@-$(CLEAN_CMD)
	@echo Clean done!

clean-all:
	@echo Cleaning all profile build directories...
	@-$(CLEAN_ALL_CMD)
	@echo Clean all done!

clean-generated:
	@echo Cleaning generated analysis outputs...
	@-$(CLEAN_GENERATED_CMD)
	@echo Clean generated done!

# Rule 'release' tự động gắn cờ LTO phù hợp với hệ điều hành
release: CXXFLAGS += $(LTO_FLAG)
release: LDFLAGS += $(LTO_FLAG)
release: clean all

native: PROFILE := native
native: clean all

portable: PROFILE := portable
portable: clean all

release-native: PROFILE := native
release-native: release

release-portable: PROFILE := portable
release-portable: release

print-config:
	@echo PROFILE=$(PROFILE)
	@echo TOOLCHAIN=$(TOOLCHAIN)
	@echo LLVM_VERSION=$(LLVM_VERSION)
	@echo LLVM_PREFIX=$(LLVM_PREFIX)
	@echo CXX=$(CXX)
	@echo TITAN_CXX23_FLOOR=$(TITAN_CXX23_FLOOR)
	@echo BUILD_DIR=$(BUILD_DIR)
	@echo ARCH_FLAGS=$(ARCH_FLAGS)
	@echo COMMON_CXX_DEFINES=$(COMMON_CXX_DEFINES)
	@echo OPENMP_FLAG=$(OPENMP_FLAG)
	@echo CXX_STDLIB_FLAG=$(CXX_STDLIB_FLAG)
	@echo CLANG_STDLIB_INCLUDE_FLAGS=$(CLANG_STDLIB_INCLUDE_FLAGS)
	@echo CLANG_OMP_INCLUDE_FLAGS=$(CLANG_OMP_INCLUDE_FLAGS)
	@echo CLANG_STDLIB_LINK_FLAGS=$(CLANG_STDLIB_LINK_FLAGS)
	@echo CLANG_STDLIB_LIBS=$(CLANG_STDLIB_LIBS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo LDFLAGS=$(LDFLAGS)

# Nạp cây phụ thuộc (Dependency Tracking)
-include $(DEPS)

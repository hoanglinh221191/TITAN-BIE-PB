# ==========================================
# TITAN PROJECT MAKEFILE (Cross-platform)
# Windows: Clang (C++23) | Linux: GCC 13 (C++23)
# Tích hợp Auto Dependency Tracking (.d)
# ==========================================

PROFILE ?= native

ifeq ($(PROFILE),native)
    ARCH_FLAGS = -march=native
else ifeq ($(PROFILE),portable)
    ARCH_FLAGS = -mavx2 -mfma -mtune=generic
else
    $(error Unsupported PROFILE '$(PROFILE)'. Use PROFILE=native or PROFILE=portable)
endif

ifeq ($(OS),Windows_NT)
    # ─── CẤU HÌNH CHO WINDOWS (CLANG) ───
    CXX  = clang++
    CC   = clang
    RM   = del /Q /F
    TARGET = titan.exe
    
    # Clang Flags: Thêm cờ chặn cảnh báo rác của Clang
    # Force Eigen to use its own aligned allocator on Windows.
    # This avoids AVX-aligned loads/stores from touching only 16-byte-aligned buffers.
    CXXFLAGS = -I. -Igmxtraj/include -O3 $(ARCH_FLAGS) -fopenmp -pipe -std=c++23 -Wall -Wno-unused-command-line-argument -MMD -MP
    CFLAGS   = -Igmxtraj/include -O3 $(ARCH_FLAGS) -ffast-math -pipe -w -Wno-incompatible-pointer-types -Wno-unused-command-line-argument -MMD -MP
    
    # Ép dùng LLD linker siêu tốc của LLVM
    LDFLAGS = -fuse-ld=lld -static -static-libgcc -static-libstdc++ -fopenmp
    
    # Clang dùng ThinLTO đa luồng
    LTO_FLAG = -flto=thin
    
    # Lệnh dọn dẹp trên Windows CMD/MSYS2
    CLEAN_CMD = $(RM) *.o *.d gmxtraj\*.o gmxtraj\*.d $(TARGET) 2>nul || exit 0
    CLEAN_GENERATED_CMD = $(RM) network_edges.csv network_nodes.csv pca_fel.csv tica_projections.csv extract_basins.sh tmp_test_dham*.in 2>nul || exit 0
else
    # ─── CẤU HÌNH CHO LINUX (GCC 13) ───
    CXX  = g++
    CC   = gcc
    RM   = rm -f
    TARGET = titan
    
    # GCC Flags: Chuẩn C++23, loại bỏ các cờ đặc thù của Clang
    # Keep Eigen on its own aligned-allocation path across Linux builds too.
    CXXFLAGS = -I. -Igmxtraj/include -O3 $(ARCH_FLAGS) -fopenmp -pipe -std=c++23 -Wall -MMD -MP
    CFLAGS   = -Igmxtraj/include -O3 $(ARCH_FLAGS) -ffast-math -pipe -w -MMD -MP
    
    # Linker mặc định của GCC (ld)
    LDFLAGS = -static-libgcc -static-libstdc++ -fopenmp
    
    # GCC dùng LTO tự động nhận diện số nhân CPU
    LTO_FLAG = -flto=auto
    
    # Lệnh dọn dẹp trên Linux bash
    CLEAN_CMD = $(RM) *.o *.d gmxtraj/*.o gmxtraj/*.d $(TARGET) 2>/dev/null || true
    CLEAN_GENERATED_CMD = $(RM) network_edges.csv network_nodes.csv pca_fel.csv tica_projections.csv extract_basins.sh tmp_test_dham*.in 2>/dev/null || true
endif

# ─── Tệp nguồn (Sources) ───
SRCS  = $(wildcard *.cpp)
OBJS  = $(SRCS:.cpp=.o)
DEPS  = $(SRCS:.cpp=.d)

CSRCS = gmxtraj/xdrfile.c gmxtraj/xdrfile_xtc.c gmxtraj/xdrfile_trr.c
COBJS = $(CSRCS:.c=.o)
CDEPS = $(CSRCS:.c=.d)

# ==========================================
# BUILD RULES
# ==========================================

.PHONY: all clean clean-generated release native portable release-native release-portable print-config

all: $(TARGET)

$(TARGET): $(OBJS) $(COBJS)
	@echo Linking $(TARGET) [$(PROFILE)]...
	@$(CXX) -o $@ $^ $(LDFLAGS)
	@echo Build successful!

%.o: %.cpp
	@echo Compiling $<...
	@$(CXX) $(CXXFLAGS) -c $< -o $@

gmxtraj/%.o: gmxtraj/%.c
	@echo Compiling $<...
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo Cleaning up...
	@-$(CLEAN_CMD)
	@echo Clean done!

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
	@echo CXX=$(CXX)
	@echo CC=$(CC)
	@echo ARCH_FLAGS=$(ARCH_FLAGS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo CFLAGS=$(CFLAGS)
	@echo LDFLAGS=$(LDFLAGS)

# Nạp cây phụ thuộc (Dependency Tracking)
-include $(DEPS) $(CDEPS)

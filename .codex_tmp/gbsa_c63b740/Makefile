# ==========================================
# TITAN PROJECT MAKEFILE (Cross-platform)
# Windows: Clang (C++23) | Linux: GCC 13 (C++23)
# Tích hợp Auto Dependency Tracking (.d)
# ==========================================

ifeq ($(OS),Windows_NT)
    # ─── CẤU HÌNH CHO WINDOWS (CLANG) ───
    CXX  = clang++
    CC   = clang
    RM   = del /Q /F
    TARGET = elec_vdw.exe
    
    # Clang Flags: Thêm cờ chặn cảnh báo rác của Clang
    CXXFLAGS = -I. -Igmxtraj/include -O3 -march=native -fopenmp -pipe -std=c++23 -Wall -Wno-unused-command-line-argument -MMD -MP
    CFLAGS   = -Igmxtraj/include -O3 -march=native -ffast-math -pipe -w -Wno-incompatible-pointer-types -Wno-unused-command-line-argument -MMD -MP
    
    # Ép dùng LLD linker siêu tốc của LLVM
    LDFLAGS = -fuse-ld=lld -static -static-libgcc -static-libstdc++ -fopenmp
    
    # Clang dùng ThinLTO đa luồng
    LTO_FLAG = -flto=thin
    
    # Lệnh dọn dẹp trên Windows CMD/MSYS2
    CLEAN_CMD = $(RM) *.o *.d gmxtraj\*.o gmxtraj\*.d $(TARGET) 2>nul || exit 0
else
    # ─── CẤU HÌNH CHO LINUX (GCC 13) ───
    CXX  = g++
    CC   = gcc
    RM   = rm -f
    TARGET = elec_vdw
    
    # GCC Flags: Chuẩn C++23, loại bỏ các cờ đặc thù của Clang
    CXXFLAGS = -I. -Igmxtraj/include -O3 -march=native -fopenmp -pipe -std=c++23 -Wall -MMD -MP
    CFLAGS   = -Igmxtraj/include -O3 -march=native -ffast-math -pipe -w -MMD -MP
    
    # Linker mặc định của GCC (ld)
    LDFLAGS = -static-libgcc -static-libstdc++ -fopenmp
    
    # GCC dùng LTO tự động nhận diện số nhân CPU
    LTO_FLAG = -flto=auto
    
    # Lệnh dọn dẹp trên Linux bash
    CLEAN_CMD = $(RM) *.o *.d gmxtraj/*.o gmxtraj/*.d $(TARGET) 2>/dev/null || true
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

all: $(TARGET)

$(TARGET): $(OBJS) $(COBJS)
	@echo Linking $(TARGET)...
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

# Rule 'release' tự động gắn cờ LTO phù hợp với hệ điều hành
release: CXXFLAGS += $(LTO_FLAG)
release: LDFLAGS += $(LTO_FLAG)
release: clean all

# Nạp cây phụ thuộc (Dependency Tracking)
-include $(DEPS) $(CDEPS)
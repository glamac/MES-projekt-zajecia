# modified "juicy make example"
# from https://makefiletutorial.com/
TARGET := messim
BUILD_DIR := ./build
SRC_DIRS := ./src
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -or -name '*.c' -or -name '*.asm')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
CPPFLAGS := $(INC_FLAGS) -MMD -MP -std=c++23 -fPIC -O3
LDFLAGS := $(shell pkg-config --static --libs glfw3) -fPIE -no-pie

$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.c.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.cpp.o : %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.asm.o: %.asm
	mkdir -p $(dir $@)
	$(CXX) -x assembler-with-cpp -c -fPIC $< -o $@

.PHONY: run
run: $(BUILD_DIR)/$(TARGET)
	./$<

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)
-include $(DEPS)

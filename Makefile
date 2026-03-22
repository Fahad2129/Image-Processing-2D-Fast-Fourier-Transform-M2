

CC       = riscv64-linux-gnu-gcc
QEMU     = qemu-riscv64
CFLAGS   = -static -O0 -g
TARGET   = fft1d

.PHONY: all run clean validate

all: $(TARGET)

$(TARGET): main.c fft_asm.s
	$(CC) $(CFLAGS) -o $@ main.c fft_asm.s -lm

run: $(TARGET)
	$(QEMU) ./$(TARGET)

validate: validate_fft.c
	gcc -O0 -o validate validate_fft.c -lm
	./validate

clean:
	rm -f $(TARGET) validate *.o *.log
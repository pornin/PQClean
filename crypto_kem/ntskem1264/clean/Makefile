# This Makefile can be used with GNU Make or BSD Make

LIB=libntskem1264_clean.a
HEADERS=additive_fft.h api.h berlekamp_massey.h bits.h ff.h m4r.h matrix_ff2.h nts_kem_errors.h nts_kem_params.h nts_kem.h polynomial.h stack.h
OBJECTS=additive_fft.o berlekamp_massey.o ff.o kem.o m4r.o matrix_ff2.o nts_kem.o polynomial.o stack.o

CFLAGS=-O3 -Wall -Wextra -Wpedantic -Werror -Wmissing-prototypes -Wredundant-decls -std=c99 -I../../../common $(EXTRAFLAGS)

all: $(LIB)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB): $(OBJECTS)
	$(AR) -r $@ $(OBJECTS)

clean:
	$(RM) $(OBJECTS)
	$(RM) $(LIB)

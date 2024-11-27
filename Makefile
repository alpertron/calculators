#
#    This file is part of Alpertron Calculators.
#
#    Copyright 2018 Dario Alejandro Alpern
#
#    Alpertron Calculators is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Alpertron Calculators is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
#
#    Use flags_coverage=--coverage in make command line to compile
#    with coverage support.
flags_factorization_1=-DFACTORIZATION_APP=1 -DFACTORIZATION_FUNCTIONS=1
flags_squares_1=
flags_other_1=
ifeq ($(bits), 32)
flags_general=-m32
else
flags_general=-D_USING64BITS_
endif
ifeq ($(flags_coverage), )
flags_general+=-g -Os
else
flags_general+=-g -O0
endif
flags_cov_and_asan=$(flags_coverage) $(flags_sanitize)
flags_factorization=$(flags_factorization_1) $(flags_cov_and_asan) $(flags_general)
flags_squares=$(flags_squares_1) $(flags_cov_and_asan) $(flags_general)
flags_other=$(flags_other_1) $(flags_cov_and_asan) $(flags_general)
h_files=batch.h bignbr.h commonstruc.h expression.h factor.h highlevel.h polynomial.h showtime.h skiptest.h
targets = ecm quad quadmod fsquares fcubes polfact dilog gaussian contfrac blockly tsqcubes sumquad divisors isprime
.PHONY : all
all: $(targets)

ecm_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c fromBlockly.c linkedbignbr.c test.c
ecm: $(ecm_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=13 $(ecm_files) -lm -o $@

blockly_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c fromBlockly.c linkedbignbr.c test.c
blockly: $(blockly_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=22 -DUSING_BLOCKLY=1 $(blockly_files) -lm -o $@

testmodmult_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c fromBlockly.c linkedbignbr.c test.c
testmultmod: $(testmodmult_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=24 $(testmodmult_files) -lm -o $@

sumquad_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c fromBlockly.c linkedbignbr.c test.c
sumquad: $(sumquad_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=25 $(sumquad_files) -lm -o $@

divisors_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c fromBlockly.c linkedbignbr.c test.c
divisors: $(divisors_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=27 $(divisors_files) -lm -o $@

gaussian_files = parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c gaussian.c GaussExpr.c bignbr.c showtime.c from_musl.c inputstr.c fft.c gcdrings.c test.c
gaussian: $(gaussian_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=12 $(gaussian_files) -lm -o $@

quad_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c quad.c quadmodLL.c bignbr.c showtime.c from_musl.c inputstr.c fft.c test.c
quad: $(quad_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=17 $(quad_files) -lm -o $@

dilog_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c dilog.c bignbr.c showtime.c from_musl.c inputstr.c fft.c test.c
dilog: $(dilog_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=11 $(dilog_files) -lm -o $@

quadmod_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
factor.c ecm.c siqs.c siqsLA.c quadmod.c quadmodLL.c bignbr.c showtime.c from_musl.c inputstr.c fft.c test.c
quadmod: $(quadmod_files) $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=16 $(quadmod_files) -lm -o $@

fsquares_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
fsquares.c tsquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c test.c
fsquares: $(fsquares_files) $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=1 $(fsquares_files) -lm -o $@

tsqcubes_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
tsqcubes.c tsquares.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c test.c
tsqcubes: $(tsqcubes_files) $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=23 $(tsqcubes_files) -lm -o $@

fcubes_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
fcubes.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c test.c
fcubes: $(fcubes_files) $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=2 $(fcubes_files) -lm -o $@

isprime_files = isprime.c test.c
isprime: $(isprime_files) $(h_files)
	gcc $(flags_other) -DDEBUG_CODE=28 $(isprime_files) -lm -o $@

contfrac_files = expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
contfrac.c bignbr.c showtime.c from_musl.c inputstr.c batch.c fft.c gcdrings.c test.c
contfrac: $(contfrac_files) $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=6 $(contfrac_files) -lm -o $@

polfact_files = expression.c parseexpr.c partition.c errors.c bigint.c linkedbignbr.c division.c baseconv.c karatsuba.c modmult.c sqroot.c \
rootseq.c lineareq.c quadraticeq.c cubiceq.c quartics.c quintics.c quinticsData.c bigrational.c output.c copyStr.c polynomial.c polyexpr.c multpoly.c divpoly.c fftpoly.c polfact.c bignbr.c showtime.c from_musl.c inputstr.c fft.c test.c intpolfact.c modpolfact.c
polfact: $(polfact_files) $(h_files)
	gcc $(flags_other) -DDEBUG_CODE=9 -DPOLYEXPR=1 $(polfact_files) -lm -o $@

clean:
	rm -f *.gcda *.gcno *.gcov $(targets)


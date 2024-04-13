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
flags_factorization_1=-D_USING64BITS_ -DFACTORIZATION_APP=1 -DFACTORIZATION_FUNCTIONS=1 -c
flags_squares_1=-D_USING64BITS_ -c
flags_other_1=-D_USING64BITS_ -c
ifeq ($(flags_coverage), )
flags_factorization=$(flags_factorization_1) -Os
flags_squares=$(flags_squares_1) -Os
flags_other=$(flags_other_1) -Os
else
flags_factorization=$(flags_factorization_1) $(flags_coverage) -g -O0
flags_squares=$(flags_squares_1) $(flags_coverage) -g -O0
flags_other=$(flags_other_1) $(flags_coverage) -g -O0
endif
h_files=batch.h bignbr.h commonstruc.h expression.h factor.h highlevel.h polynomial.h showtime.h skiptest.h
targets = ecm quad quadmod fsquares fcubes polfact dilog contfrac blockly tsqcubes sumquad divisors isprime
.PHONY : all
all: $(targets)

ecm: expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco ecmfront.fco sumSquares.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco fromBlockly.fco linkedbignbr.fco test_ecm.fco
	gcc expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco ecmfront.fco sumSquares.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco fromBlockly.fco linkedbignbr.fco test_ecm.fco $(flags_coverage) -lm -o $@

blockly: expression.fbo parseexpr.fbo partition.fbo errors.fbo copyStr.fbo bigint.fbo division.fbo baseconv.fbo karatsuba.fbo modmult.fbo sqroot.fbo \
factor.fbo ecm.fbo siqs.fbo siqsLA.fbo ecmfront.fbo sumSquares.fco bignbr.fbo showtime.fbo from_musl.fbo inputstr.fbo batch.fbo fft.fbo gcdrings.fbo fromBlockly.fbo linkedbignbr.fbo test_blockly.fbo
	gcc expression.fbo parseexpr.fbo partition.fbo errors.fbo copyStr.fbo bigint.fbo division.fbo baseconv.fbo karatsuba.fbo modmult.fbo sqroot.fbo \
factor.fbo ecm.fbo siqs.fbo siqsLA.fbo ecmfront.fbo sumSquares.fco bignbr.fbo showtime.fbo from_musl.fbo inputstr.fbo batch.fbo fft.fbo gcdrings.fbo fromBlockly.fbo linkedbignbr.fbo test_blockly.fbo $(flags_coverage) -lm -o $@

sumquad: expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco ecmfront.fco sumSquares.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco fromBlockly.fco linkedbignbr.fco test_sumquad.fco
	gcc expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco ecmfront.fco sumSquares.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco fromBlockly.fco linkedbignbr.fco test_sumquad.fco $(flags_coverage) -lm -o $@

divisors: expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco ecmfront.fco sumSquares.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco fromBlockly.fco linkedbignbr.fco test_divisors.fco
	gcc expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco ecmfront.fco sumSquares.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco fromBlockly.fco linkedbignbr.fco test_divisors.fco $(flags_coverage) -lm -o $@

gaussian: parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco gaussian.fco GaussExpr.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco gcdrings.fco test_gaussian.fco
	gcc parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco gaussian.fco GaussExpr.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco gcdrings.fco test_gaussian.fco $(flags_coverage) -lm -o $@

quad: expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco quad.fco quadmodLL.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quad.fco
	gcc expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco quad.fco quadmodLL.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quad.fco $(flags_coverage) -lm -o $@

dilog: expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco dilog.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_dilog.fco
	gcc expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco dilog.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_dilog.fco $(flags_coverage) -lm -o $@

quadmod: expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco quadmod.fco quadmodLL.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quadmod.fco
	gcc expression.fco parseexpr.fco partition.fco errors.fco copyStr.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco siqsLA.fco quadmod.fco quadmodLL.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quadmod.fco $(flags_coverage) -lm -o $@

fsquares: expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo tsquares.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fsquares.sqo
	gcc expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo tsquares.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fsquares.sqo $(flags_coverage) -lm -o $@

tsqcubes: expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
tsqcubes.sqo tsquares.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_tsqcubes.sqo
	gcc expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
tsqcubes.sqo tsquares.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_tsqcubes.sqo $(flags_coverage) -lm -o $@

fcubes: expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fcubes.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fcubes.sqo
	gcc expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fcubes.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fcubes.sqo $(flags_coverage) -lm -o $@

isprime: isprime.sqo test_isprime.sqo
	gcc isprime.sqo test_isprime.sqo $(flags_coverage) -lm -o $@

contfrac: expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_contfrac.sqo
	gcc expression.sqo parseexpr.sqo partition.sqo errors.sqo copyStr.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_contfrac.sqo $(flags_coverage) -lm -o $@

polfact: expression.oto parseexpr.plo partition.oto errors.oto bigint.oto linkedbignbr.oto division.oto baseconv.oto karatsuba.oto modmult.oto sqroot.oto \
rootseq.oto lineareq.oto quadraticeq.oto cubiceq.oto quartics.oto quintics.oto quinticsData.oto bigrational.oto output.oto copyStr.oto polynomial.oto polyexpr.oto multpoly.oto divpoly.oto fftpoly.oto polfact.oto bignbr.oto showtime.oto from_musl.oto inputstr.oto fft.oto test_polfact.oto intpolfact.oto modpolfact.oto
	gcc expression.oto parseexpr.plo partition.oto errors.oto bigint.oto linkedbignbr.oto division.oto baseconv.oto karatsuba.oto modmult.oto sqroot.oto \
rootseq.oto lineareq.oto quadraticeq.oto cubiceq.oto quartics.oto quintics.oto quinticsData.oto bigrational.oto output.oto copyStr.oto polynomial.oto polyexpr.oto multpoly.oto divpoly.oto fftpoly.oto polfact.oto bignbr.oto showtime.oto from_musl.oto inputstr.oto fft.oto test_polfact.oto intpolfact.oto modpolfact.oto $(flags_coverage) -lm -o $@

%.fco : %.c $(h_files)
	gcc $(flags_factorization) -o $@ $<

%.fbo : %.c $(h_files)
	gcc $(flags_factorization) -DUSING_BLOCKLY=1 -o $@ $<

%.oto : %.c $(h_files)
	gcc $(flags_other) -o $@ $<

%.plo : %.c $(h_files)
	gcc -DPOLYEXPR=1 $(flags_other) -o $@ $<

%.sqo : %.c $(h_files)
	gcc $(flags_squares) -o $@ $<

test_ecm.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=13 -o $@ $<
	
test_blockly.fbo: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=22 -DUSING_BLOCKLY=1 -o $@ $<
	
test_sumquad.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=25 -o $@ $<
  
test_divisors.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=27 -o $@ $<
  
test_gaussian.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=12 -o $@ $<
	
test_quadmod.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=16 -o $@ $<
	
test_quad.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=17 -o $@ $<
	
test_dilog.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=11 -o $@ $<
	
test_fsquares.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=1 -o $@ $<
	
test_tsqcubes.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=23 -o $@ $<
	
test_fcubes.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=2 -o $@ $<
	
test_contfrac.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=6 -o $@ $<
	
test_polfact.oto: test.c $(h_files)
	gcc $(flags_other) -DDEBUG_CODE=9 -o $@ $<

test_isprime.sqo: test.c $(h_files)
	gcc $(flags_other) -DDEBUG_CODE=28 -o $@ $<

clean:
	rm -f *.oto *.fco *.fbo *.sqo *.plo *.gcda *.gcno *.gcov $(targets)


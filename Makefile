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
flags_factorization=-D_USING64BITS_ -DFACTORIZATION_APP=1 -c -Os
flags_squares=-D_USING64BITS_ -DFSQUARES_APP=1 -c -Os
flags_other=-D_USING64BITS_ -c -Os
h_files=batch.h bignbr.h commonstruc.h expression.h factor.h highlevel.h polynomial.h showtime.h skiptest.h
targets = ecm quad quadmod fsquares fcubes polfact dilog contfrac
.PHONY : all
all: $(targets)

ecm: expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco ecmfront.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco test_ecm.fco
	gcc expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco ecmfront.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco batch.fco fft.fco gcdrings.fco test_ecm.fco -lm -o $@

quad: expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco quad.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quad.fco
	gcc expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco quad.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quad.fco -lm -o $@

dilog: expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco dilog.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_dilog.fco
	gcc expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco dilog.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_dilog.fco -lm -o $@

quadmod: expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco quadmod.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quadmod.fco
	gcc expression.fco partition.fco errors.fco bigint.fco division.fco baseconv.fco karatsuba.fco modmult.fco sqroot.fco \
factor.fco ecm.fco siqs.fco quadmod.fco bignbr.fco showtime.fco from_musl.fco inputstr.fco fft.fco test_quadmod.fco -lm -o $@

fsquares: expression.sqo partition.sqo errors.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo fcubes.sqo contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fsquares.sqo
	gcc expression.sqo partition.sqo errors.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo fcubes.sqo contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fsquares.sqo -lm -o $@

fcubes: expression.sqo partition.sqo errors.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo fcubes.sqo contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fcubes.sqo
	gcc expression.sqo partition.sqo errors.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo fcubes.sqo contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_fcubes.sqo -lm -o $@

contfrac: expression.sqo partition.sqo errors.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo fcubes.sqo contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_contfrac.sqo
	gcc expression.sqo partition.sqo errors.sqo bigint.sqo division.sqo baseconv.sqo karatsuba.sqo modmult.sqo sqroot.sqo \
fsquares.sqo fcubes.sqo contfrac.sqo bignbr.sqo showtime.sqo from_musl.sqo inputstr.sqo batch.sqo fft.sqo gcdrings.sqo test_contfrac.sqo -lm -o $@

polfact: expression.oto partition.oto errors.oto bigint.oto linkedbignbr.oto division.oto baseconv.oto karatsuba.oto modmult.oto sqroot.oto \
rootseq.oto quintics.oto bigrational.oto output.oto polynomial.oto polyexpr.oto multpoly.oto divpoly.oto fftpoly.oto polfact.oto bignbr.oto showtime.oto from_musl.oto inputstr.oto fft.oto test_polfact.oto intpolfact.oto
	gcc expression.oto partition.oto errors.oto bigint.oto linkedbignbr.oto division.oto baseconv.oto karatsuba.oto modmult.oto sqroot.oto \
rootseq.oto quintics.oto bigrational.oto output.oto polynomial.oto polyexpr.oto multpoly.oto divpoly.oto fftpoly.oto polfact.oto bignbr.oto showtime.oto from_musl.oto inputstr.oto fft.oto test_polfact.oto intpolfact.oto -lm -o $@

%.fco : %.c $(h_files)
	gcc $(flags_factorization) -o $@ $<

%.oto : %.c $(h_files)
	gcc $(flags_other) -o $@ $<

%.sqo : %.c $(h_files)
	gcc $(flags_squares) -o $@ $<

test_ecm.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=13 -o $@ $<
	
test_quadmod.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=16 -o $@ $<
	
test_quad.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=17 -o $@ $<
	
test_dilog.fco: test.c $(h_files)
	gcc $(flags_factorization) -DDEBUG_CODE=11 -o $@ $<
	
test_fsquares.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=1 -o $@ $<
	
test_fcubes.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=2 -o $@ $<
	
test_contfrac.sqo: test.c $(h_files)
	gcc $(flags_squares) -DDEBUG_CODE=6 -o $@ $<
	
test_polfact.oto: test.c $(h_files)
	gcc $(flags_other) -DDEBUG_CODE=9 -o $@ $<

clean:
	rm -f *.oto *.fco *.sqo ecm $(targets)


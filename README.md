I have written 14 calculators in my free time. Their source code is open source, so you can use the code for your own projects.

I started writing these calculators as Java applets in 1997. In 2015 I ported them to C language. The code is compiled to WebAssembly or asm.js using Emscripten so it can run inside Web browsers.

The complete source code size is about 70 000 lines of code.

### Running the programs as standalone executables

The main purpose of compiling the code is to run coverage tests. The file ``coverage.out.old`` holds the results of running all calculators with different inputs.
After making a change in the code, I run this test again and the output should not be changed. When adding a new feature to the calculators, I add the new case to the coverage, and change the results file if the results are correct.
These tests use the files ``Makefile`` and ``coverage``.

You can use ``Makefile`` to generate standalone executables. Just run ``make clean`` and then ``make``.
If you want to build only one of the calculators, you can run ``make calculator``, where ``calculator`` is one of the following words:

-   ``ecm`` (integer factorization)
-   ``gaussian`` (Gaussian integer factorization)
-   ``quad`` (two-variable quadratic integer solver)
-   ``dilog`` (discrete logarithm solver)
-   ``quadmod`` (quadratic modular equation solver)
-   ``fsquares`` (decompose number in sum of squares)
-   ``fcubes`` (decompose number in sum of cubes)
-   ``tsqcubes`` (decompose number in sum of two squares and a cube, fifth or seventh power)
-   ``contfrac`` (continued fraction calculator)
-   ``polfact`` (polynomial equation solver and factorization)
-   ``sumquad`` (all sums of two squares)

After building the calculators, you can run them. For example:

``./polfact 0 "x^2+x+3" 2``

The output is in HTML and it starts with a digit (this is used by the JavaScript code).
All calculators can process expressions.

If you run the calculator without command line parameters, the program will show the expected arguments.

### Static code analysis

Software quality assurance is realized with the coverage test explained in the previous section and static analysis. The latter is done by tools that scan the source code. The programs do not run in this case.

Results of static analysis and code coverage of this software using [Codacy](https://app.codacy.com/gh/alpertron/calculators/dashboard).

Badges from [SonarCloud](https://sonarcloud.io/summary/overall?id=alpertron_calculators):

[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=bugs)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Code Smells](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=code_smells)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=coverage)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Duplicated Lines (%)](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=duplicated_lines_density)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=ncloc)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Maintainability Rating](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=sqale_rating)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=alert_status)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Reliability Rating](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=reliability_rating)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Security Rating](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=security_rating)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Technical Debt](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=sqale_index)](https://sonarcloud.io/dashboard?id=alpertron_calculators)
[![Vulnerabilities](https://sonarcloud.io/api/project_badges/measure?project=alpertron_calculators&metric=vulnerabilities)](https://sonarcloud.io/dashboard?id=alpertron_calculators)

If you like these calculators and you want to support free software, you can donate via Paypal by clicking in the button below:

[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=MR65QPWZM5JT6)

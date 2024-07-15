.PHONY: all s21_matrix_oop.a test gcov_report check leaks valgrind rebuild clean clean_test

CC = g++ -m64 -Wall -Wextra -Werror -std=c++17
LIBS = -I/usr/local/include/ -L/usr/local/lib -lgtest -lgmock
HEADER = s21_matrix_oop.h
SRC = s21_matrix.cpp
TESTS = tests.cpp

all: clean s21_matrix_oop.a test gcov_report

s21_matrix_oop.a: $(SRC) $(HEADER)
	$(CC) -c $(SRC) $(HEADER)
	ar rcs s21_matrix_oop.a *.o
	ranlib s21_matrix_oop.a
	rm -rf *.o

test: s21_matrix_oop.a tests.cpp
	$(CC) $(TESTS) s21_matrix_oop.a -o test.out $(LIBS)
	./test.out

gcov_report: clean_test
	$(CC) --coverage $(TESTS) $(SRC) $(HEADER) -lstdc++ -o test.out $(LIBS)
	./test.out
	lcov -t "./test.out" --ignore-errors source -o report.info --no-external -c -d .
	genhtml --ignore-errors source -o report report.info
	open ./report/index.html

check:
	cppcheck --enable=all --suppress=missingIncludeSystem ./
	clang-format -i *.cpp
	clang-format -i *.h

leaks:
	$(CC) $(SRC) $(HEADER) tests.cpp -o main.out $(LIBS)
	leaks --atExit -- ./main.out
	rm -rf main.out

valgrind:
	$(CC) $(SRC) $(HEADER) tests.cpp -o main.out $(LIBS)
	valgrind --tool=memcheck --leak-check=yes ./main.out
	rm -f main.out

rebuild: clean all

clean: clean_test
	rm -rf ./report

clean_test:
	rm -f *.o *.a *.gcno *.gcda *.info *.out *.gcov *.gch

#! /bin/bash
#### this script will find all the test_class modules in testing/ directory
#### and generate the unit testing executable folly
#### folly usage is as follows:
#### ./folly (will run all unit tests)
#### ./folly class_name (will run all unit test of a particular class)
#### ./folly class_name unit_test (will run a particular unit test for class) 
#### ./folly list (will list all testable class_names)
#### ./folly class_name list (will list all unit tests for particular class)


if [ "$#" -eq 0 ]; then
# run all unit tests
echo "run all unit test - place holder"
fi

if [ "$#" -eq 1 ]; then
    if [ "$1" == "list" ]; then
echo "List of testable class names."
echo "-----------------------------"
# return all class names in testing/ directory with filename test_classname.f90
	ls testing/ | grep "^test_" | sed "s/test_//g"| sed "s/.f90//g"
echo "-----------------------------"

    fi
fi
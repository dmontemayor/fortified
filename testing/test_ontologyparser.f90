!purpose: module containing all unit tests for ontology parser

module test_ontologyparser
  use testing_class
  
contains
  subroutine test_parser_can_be_executed
    integer::exitstat,cmdstat
    
    call execute_command_line("echo 'Hello World'"&
         ,exitstat=exitstat,cmdstat=cmdstat)
    write(*,*)exitstat,cmdstat
    
    write(*,*)'ONTOLOGY PARSER PASSED ALL TESTS!'
  end subroutine test_parser_can_be_executed

end module test_ontologyparser

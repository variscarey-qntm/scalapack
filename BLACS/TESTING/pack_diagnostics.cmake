execute_process(
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
          ${TEST_PROG} ${MODE} ${MPIEXEC_POSTFLAGS}
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE diagnostic
  TIMEOUT 30)

if("${result}" STREQUAL "0" OR "${result}" MATCHES "[Tt]imeout")
  message(FATAL_ERROR "Expected MPI_Abort, got ${result}: ${diagnostic}")
endif()
if("${diagnostic}" MATCHES "TEST FAILURE")
  message(FATAL_ERROR "${diagnostic}")
endif()
foreach(required
    "BLACS ERROR: ${OPERATION} failed"
    "MPI error [0-9]+"
    "Contxt=[0-9]+, pnum=0, grid=1x1, coordinates=\\{0,0\\}"
    "scope=a, scope_rank=0, scope_size=1"
    "count=1, datatype_handle=0x"
    "on line [0-9]+ of file")
  if(NOT "${diagnostic}" MATCHES "${required}")
    message(FATAL_ERROR "Missing '${required}': ${diagnostic}")
  endif()
endforeach()
if("${MODE}" STREQUAL "size")
  if(NOT "${diagnostic}" MATCHES "packed_size=unavailable, buffer=not allocated")
    message(FATAL_ERROR "Missing sizing diagnostics: ${diagnostic}")
  endif()
else()
  if(NOT "${diagnostic}" MATCHES
      "input=.*, buffer=.*, buffer_capacity=[0-9]+ bytes, outstanding_operations=0, initial_position=0")
    message(FATAL_ERROR "Missing packing diagnostics: ${diagnostic}")
  endif()
endif()
if("${MODE}" STREQUAL "error-string" AND
   NOT "${diagnostic}" MATCHES "MPI error text unavailable")
  message(FATAL_ERROR "Missing fallback error text: ${diagnostic}")
endif()

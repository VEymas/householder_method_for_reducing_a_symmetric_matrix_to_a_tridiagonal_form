target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE
  householder_for_sym_matrix.cpp
)

if(BUILD_BLAS)
  find_package(BLAS REQUIRED)
  target_link_libraries(${CMAKE_PROJECT_NAME}
    ${BLAS_LIBRARIES}
  )
  target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE USE_BLAS)
else()
  target_link_libraries(${CMAKE_PROJECT_NAME})
endif()

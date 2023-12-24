add_executable(main
  main.cpp
)

target_link_libraries(main
${CMAKE_PROJECT_NAME}
)
add_executable(fpotencia_test
	test-models.cpp

	nr-current-solver-test.cpp
	nr-polar-solver-test.cpp
)

target_link_libraries(fpotencia_test
	PUBLIC
		Catch2::Catch2WithMain
		fpotencia
)

include(Catch)
catch_discover_tests(fpotencia_test)

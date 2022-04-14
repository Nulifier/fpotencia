add_executable(fpotencia_test
	test.cpp
	test-models.cpp
)

target_link_libraries(fpotencia_test
	PUBLIC
		fpotencia
		Catch2::Catch2WithMain
)

include(CTest)
include(Catch)

catch_discover_tests(fpotencia_test)

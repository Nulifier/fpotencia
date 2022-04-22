add_executable(fpotencia_profile
	../test/test-models.cpp

	profile.cpp
)

target_include_directories(fpotencia_profile
	PRIVATE
		${CMAKE_SOURCE_DIR}/test
)

target_link_libraries(fpotencia_profile
	PRIVATE
		fpotencia
)

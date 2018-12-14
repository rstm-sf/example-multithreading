# The function adds the target of the executable file.
function(custom_add_executable_from_dir TARGET)
    file(GLOB TARGET_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    )
    add_executable(${TARGET} ${TARGET_SRC})
    target_link_libraries(${TARGET}
        GTest::GTest
        GTest::Main
    )
endfunction()
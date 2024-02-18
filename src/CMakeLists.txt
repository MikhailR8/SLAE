#Записываем в переменную SOURCES list всех .cpp и .hpp из директории и субдиректорий
file(GLOB_RECURSE SOURCES *.hpp *.cpp)
add_library(${PROJECT_NAME} STATIC ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
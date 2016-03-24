set(QYates_SRCS QYates.cpp main.cpp)
set(QYates_H QYates.h)
set(QYates_UIS QYates.ui)
set(QYates_RCC icons.qrc)

# start building interface from moc
QT4_WRAP_CPP(QYates_MOC_H ${QYates_H})
QT4_WRAP_UI(QYates_UIS_H ${QYates_UIS})
QT4_ADD_RESOURCES(QYates_RCC_SRCS ${QYates_RCC})


source_group("Resources" FILES
  ${QYates_UIS}
)
 
source_group("Generated" FILES
  ${QYates_UIS_H}
  ${QYates_MOC_H}
)

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS} --debug)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(QYates ${QYates_SRCS} ${QYates_MOC_H} ${QYates_UIS_H} ${QYates_RCC_SRCS})

if(WIN32)
  target_link_libraries(QYates ${SCIENTIFIC_LIB_DIR}/libscientific.dll.a ${QT_LIBRARIES} "-mwindows") # "-mwindows" serve per non aprire la cmd di windows.
else()
  target_link_libraries(QYates m scientific ${QT_LIBRARIES})
endif()

install(TARGETS QYates DESTINATION bin COMPONENT ApplicationData)



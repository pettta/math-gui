include(ExternalProject)

if(TARGET boost_ep)
  # Boost external project already configured.
  return()
endif()

set(BOOST_EXTERNAL_VERSION 1.86.0)
string(REPLACE "." "_" BOOST_EXTERNAL_VERSION_UNDERSCORE "${BOOST_EXTERNAL_VERSION}")
set(BOOST_EXTERNAL_BASE_URL "https://archives.boost.io/release")
set(BOOST_EXTERNAL_ARCHIVE "${BOOST_EXTERNAL_BASE_URL}/${BOOST_EXTERNAL_VERSION}/source/boost_${BOOST_EXTERNAL_VERSION_UNDERSCORE}.tar.gz")

set(BOOST_EXTERNAL_PREFIX "${CMAKE_BINARY_DIR}/third_party/boost")
set(BOOST_EXTERNAL_INSTALL_DIR "${BOOST_EXTERNAL_PREFIX}/install")

ExternalProject_Add(boost_ep
  PREFIX "${BOOST_EXTERNAL_PREFIX}"
  URL "${BOOST_EXTERNAL_ARCHIVE}"
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  UPDATE_DISCONNECTED YES
  BUILD_IN_SOURCE YES
  CONFIGURE_COMMAND ./bootstrap.sh --prefix=<INSTALL_DIR> --with-libraries=filesystem,system,random
  BUILD_COMMAND ./b2 install --prefix=<INSTALL_DIR> --layout=system variant=release link=static threading=multi cxxflags=-std=c++17
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(boost_ep install_dir)
set(BOOST_EXTERNAL_INSTALL_DIR "${install_dir}" CACHE PATH "Boost external install directory" FORCE)
set(BOOST_EXTERNAL_INCLUDE_DIR "${BOOST_EXTERNAL_INSTALL_DIR}/include" CACHE PATH "Boost external include directory" FORCE)
set(BOOST_EXTERNAL_LIBRARY_DIR "${BOOST_EXTERNAL_INSTALL_DIR}/lib" CACHE PATH "Boost external library directory" FORCE)

file(MAKE_DIRECTORY "${BOOST_EXTERNAL_INCLUDE_DIR}")
file(MAKE_DIRECTORY "${BOOST_EXTERNAL_LIBRARY_DIR}")

add_library(Boost::headers INTERFACE IMPORTED)
set_target_properties(Boost::headers PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${BOOST_EXTERNAL_INCLUDE_DIR}")
add_dependencies(Boost::headers boost_ep)

foreach(_boost_lib IN ITEMS system filesystem random)
  add_library(Boost::${_boost_lib} STATIC IMPORTED)
  set_target_properties(Boost::${_boost_lib} PROPERTIES
    IMPORTED_LOCATION "${BOOST_EXTERNAL_LIBRARY_DIR}/libboost_${_boost_lib}.a"
    INTERFACE_INCLUDE_DIRECTORIES "${BOOST_EXTERNAL_INCLUDE_DIR}"
  )
  target_link_libraries(Boost::${_boost_lib} INTERFACE Boost::headers)
  add_dependencies(Boost::${_boost_lib} boost_ep)
endforeach()

target_link_libraries(Boost::filesystem INTERFACE Boost::system)

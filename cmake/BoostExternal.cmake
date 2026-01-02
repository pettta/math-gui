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
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(boost_ep source_dir)
set(BOOST_EXTERNAL_INSTALL_DIR "${source_dir}" CACHE PATH "Boost external install directory" FORCE)
set(BOOST_EXTERNAL_INCLUDE_DIR "${source_dir}" CACHE PATH "Boost external include directory" FORCE)

add_library(Boost::headers INTERFACE IMPORTED)
set_target_properties(Boost::headers PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${BOOST_EXTERNAL_INCLUDE_DIR}")
add_dependencies(Boost::headers boost_ep)

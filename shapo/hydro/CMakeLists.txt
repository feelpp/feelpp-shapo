###  CMakeLists.txt; coding: utf-8 ---

#  Author(s):<christophe.prudhomme@feelpp.org>
#       Date: 06 Feb 2015
#
#  Copyright (C) 2015-2017 Feel++ Consortium
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
cmake_minimum_required(VERSION 3.21)

add_library(feelpp_shapo_hydro SHARED stokes.cpp)
target_compile_definitions(feelpp_shapo_hydro PRIVATE FEELPP_DIM=3)
target_include_directories(feelpp_shapo_hydro PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
  $<INSTALL_INTERFACE:include/feelpp> ) 
target_link_libraries(feelpp_shapo_hydro PUBLIC Feelpp::feelpp)

feelpp_add_application( hydro_stokes SRCS stokes.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=3 LINK_LIBRARIES feelpp_shapo_hydro EXEC STOKES_SHAPE_OPT INSTALL  )
feelpp_add_testcase( cases PREFIX shapo CATEGORY hydro )

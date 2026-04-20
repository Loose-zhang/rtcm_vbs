cmake_minimum_required(VERSION 3.10)
project(rtcm_vbs_hp C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT WIN32)
    add_compile_options(-Wall -Wno-unused-parameter -Wno-unused-function
                        -Wno-unused-variable -Wno-unused-but-set-variable)
endif()

set(RTKLIB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(RTKLIB_SRC  ${RTKLIB_ROOT}/src)
set(RTKLIB_INC  ${RTKLIB_ROOT}/include)

# Minimal RTKLIB sources we actually need for RTCM VBS. Keeps compile fast
# and avoids pulling in GUI/stream/receiver-decoder code.
set(RTKLIB_SOURCES
    ${RTKLIB_SRC}/rtcm.c
    ${RTKLIB_SRC}/rtcm2.c
    ${RTKLIB_SRC}/rtcm3.c
    ${RTKLIB_SRC}/rtcm3e.c
    ${RTKLIB_SRC}/ephemeris.c
    ${RTKLIB_SRC}/rtkcmn.c
    ${RTKLIB_SRC}/preceph.c
    ${RTKLIB_SRC}/sbas.c
    ${RTKLIB_SRC}/options.c
    ${RTKLIB_SRC}/trace.c
    ${RTKLIB_SRC}/tides.c
    ${RTKLIB_SRC}/ionex.c
    ${RTKLIB_SRC}/datum.c
    ${RTKLIB_SRC}/geoid.c
    ${RTKLIB_SRC}/B2b.c
    ${RTKLIB_SRC}/f2c/cal2jd.c
    ${RTKLIB_SRC}/f2c/dat.c
    ${RTKLIB_SRC}/f2c/dehanttideinel.c
    ${RTKLIB_SRC}/f2c/gmf.c
    ${RTKLIB_SRC}/f2c/gpt.c
    ${RTKLIB_SRC}/f2c/norm8.c
    ${RTKLIB_SRC}/f2c/sprod.c
    ${RTKLIB_SRC}/f2c/st1idiu.c
    ${RTKLIB_SRC}/f2c/st1isem.c
    ${RTKLIB_SRC}/f2c/st1l1.c
    ${RTKLIB_SRC}/f2c/step2diu.c
    ${RTKLIB_SRC}/f2c/step2lon.c
    ${RTKLIB_SRC}/f2c/vmf1.c
    ${RTKLIB_SRC}/f2c/vmf1_ht.c
    ${RTKLIB_SRC}/f2c/zero_vec8.c
)

add_library(rtklib_min STATIC ${RTKLIB_SOURCES})
target_include_directories(rtklib_min PUBLIC ${RTKLIB_INC})
target_compile_definitions(rtklib_min PUBLIC
    ENAGLO ENAGAL ENACMP ENAQZS
)
# RTKLIB sources occasionally rely on implicit declarations of stdio funcs
# and contain && / || parenthesisation that newer compilers warn about.
# Silence those to let upstream code compile unchanged.
if(NOT MSVC)
    target_compile_options(rtklib_min PRIVATE
        -Wno-implicit-function-declaration
        -Wno-logical-op-parentheses
        -Wno-parentheses
        -Wno-format
        -Wno-incompatible-pointer-types
        -Wno-int-conversion
        -Wno-deprecated-non-prototype
        -Wno-pointer-sign
    )
endif()
if(UNIX)
    target_link_libraries(rtklib_min PUBLIC m pthread)
elseif(WIN32)
    target_link_libraries(rtklib_min PUBLIC ws2_32 winmm)
endif()

add_executable(rtcm_vbs_hp
    main.c
    vbs_core.c
    vbs_merge.c
    tcp_io.c
    net_interp.c
    rtklib_glue.c
)

target_link_libraries(rtcm_vbs_hp PRIVATE rtklib_min)

# run() holds several channel_t, each with a full rtcm_t (nav_t is huge). The
# default 1 MB stack on Windows MSVC is too small and the process exits instantly.
# MinGW/Clang on Windows uses the same default; match MSVC stack size.
if(MSVC)
    set_target_properties(rtcm_vbs_hp PROPERTIES LINK_FLAGS "/STACK:8388608")
elseif(WIN32 AND NOT MSVC)
    set_target_properties(rtcm_vbs_hp PROPERTIES LINK_FLAGS "-Wl,--stack,8388608")
endif()

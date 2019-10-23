workspace "NN"
    architecture "x64"
    configurations
    {
        "Debug",
        "Release",
    }

    staticruntime "on"
    language "C"
    objdir "obj/%{cfg.buildcfg}"
    targetdir "bin/%{cfg.buildcfg}"
    debugdir "build"

project "NN"
    kind "StaticLib"

    files
    {
        "delaunay.c",
        "hash.c",
        "istack.c",
        "lpi.c",
        "minell.c",
        "nnai.c",
        "nnpi.c",
        "nncommon.c",
        "nncommon-vulnerable.c",
        "triangle.c",
    }
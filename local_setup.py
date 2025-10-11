import subprocess
import sys
import os 

print(f"Running local setup on {sys.platform}")


def run_cmake_build():
    subprocess.run(["cmake", "-B", "build", "."])
    subprocess.run(["cmake", "--build", "build", "--config", "Release"])


if sys.platform.startswith("darwin"):
    try:
        run_cmake_build()
        os.chdir("bin")
        print("Now you can run the executables! Running Metal Engine...")
        subprocess.run(["./metal_engine"])

    except subprocess.CalledProcessError as e:
        print(f"Error: cmake exited with {e.returncode}")
        sys.exit(e.returncode)


if sys.platform.startswith("linux"):
    try:
        run_cmake_build()
        os.chdir("bin")
        print("Now you can run the executables! Running Vulkan Engine...")
        subprocess.run(["./vulkan_engine"])

    except subprocess.CalledProcessError as e:
        print(f"Error: cmake exited with {e.returncode}")
        sys.exit(e.returncode)

        
if sys.platform.startswith("win"):
    # Need the path prefix to account for MSVC path installation
    # To find this manually run "echo %PATH%" in a Visual Studio Developer Command Prompt & look at the start and end 
    PATH_PREFIX = r"C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\VC\Linux\bin\ConnectionManagerExe;C:\Program Files\Microsoft Visual Studio\2022\Community\VC\vcpkg;C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\bin\HostX64\x64;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\VC\VCPackages;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\TestWindow;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\TeamFoundation\Team Explorer;C:\Program Files\Microsoft Visual Studio\2022\Community\MSBuild\Current\bin\Roslyn;C:\Program Files\Microsoft Visual Studio\2022\Community\Team Tools\DiagnosticsHub\Collector;C:\Program Files (x86)\Windows Kits\10\bin\10.0.26100.0\\x64;C:\Program Files (x86)\Windows Kits\10\bin\\x64;C:\Program Files\Microsoft Visual Studio\2022\Community\\MSBuild\Current\Bin\amd64;C:\Windows\Microsoft.NET\Framework64\v4.0.30319;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\;C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\;"
    print(f"Running on Windows w/ Path={PATH_PREFIX + os.pathsep + os.environ['PATH']}")
    os.environ["PATH"] = PATH_PREFIX + os.pathsep + os.environ["PATH"]
    try:
        subprocess.run(["cmake", "-B", "build", "."])
        subprocess.run(["cmake", "--build", "build", "--config", "Release"])
        os.chdir(os.path.join("bin", "Release"))
        print("Now you can run the executables! Running Vulkan Engine...")
        subprocess.run(["vulkan_engine.exe"])

    except subprocess.CalledProcessError as e:
        print(f"Error: cmake exited with {e.returncode}")
        sys.exit(e.returncode)


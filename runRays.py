#!/usr/bin/python3

if __name__ == "__main__":
    import sys
    sys.path.insert(0, "RayTrace")
    sys.path.insert(0, "Atmosphere")
    sys.path.append("./Env")
    sys.path.insert(0, "input")
    sys.path.insert(0, "venv/lib/python3.7/site-packages")
    print(sys.path)
    import RayTrace
    RayTrace.main()

# Wi-Fi Scheduling Solution

1. The scenario is defined via links; Extend `Class LinkBase` to define the *average link rate*, `ThruApp/RTApp/DLApp` for each AC groups;
2. The problem is defnied in `demo.py` file, which iterate over the provided links definition and formulate a CVXPY problem;

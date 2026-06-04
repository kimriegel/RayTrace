from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import sys
sys.path.insert(0, "Env/")
sys.path.insert(0, "input")
sys.path.insert(0, "venv/lib/python3.7/site-packages")
sys.path.insert(0, "RayTrace/")
print(sys.path)
from newparameterfile import create_parameters

app = FastAPI()



app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # allows everything (fine for testing)
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SimulationParameters(BaseModel):
    Fs: float
    xinitial: float
    yinitial: float
    boomspacing: float
    h: float


@app.post("/run-simulation")
def run_simulation(params: SimulationParameters):
    simulation_params = create_parameters(
        Fs=params.Fs,
        xinitial=params.xinitial,
        yinitial=params.yinitial,
        boomspacing=params.boomspacing,
        h=params.h,
    )
    print(simulation_params)
    import RayTrace
    RayTrace.main()

    return {
        "message": "Simulation parameters were successfully sent to Python.",
        "parameters": simulation_params
    }



# #!/usr/bin/python3

# if __name__ == "__main__":
#     import sys
#     sys.path.insert(0, "RayTrace")

#     import RayTrace
#     RayTrace.main()

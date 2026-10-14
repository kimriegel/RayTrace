from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import sys


sys.path.insert(0, "Env/")
sys.path.insert(0, "input")
sys.path.insert(0, "venv/lib/python3.7/site-packages")
sys.path.insert(0, "RayTrace/")
print(sys.path)

from newparameterfile import set_parameters

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
    zinitial: float
    radius: float
    soundspeed: float
    ps: float
    Temp: float
    hr: float
    theta: float
    phi: float
    boomspacing: float
    h: float


@app.post("/run-simulation")
def run_simulation(params: SimulationParameters):

    print("Values received from website:")
    print("Fs:", params.Fs)
    print("xinitial:", params.xinitial)
    print("boomspacing:", params.boomspacing)
    print("h:", params.h)

    set_parameters(
        params.Fs,
        params.xinitial,
        params.yinitial,
        params.zinitial,
        params.radius,
        params.soundspeed,
        params.ps,
        params.Temp,
        params.hr,
        params.theta,
        params.phi,
        params.boomspacing,
        params.h
    )
    

    import RayTrace

    print(RayTrace)
    print(RayTrace.__file__)
    print(dir(RayTrace))

    graph_data = RayTrace.main()   
    #Run simulation, when RayTrace gives graph information back, save it in variable graph_data

    # Check whether the user stopped the simulation
    if graph_data is None:
        return {
            "message": "Simulation stopped by user.",
            "parameters": params,
            "graph_data": None,
            "stopped": True
        }

    return {
        "message": "Simulation parameters were successfully ran and sent to Python.",
        "parameters": params,
        "graph_data": graph_data,
        "stopped": False   #FIXXXXXXXXX
}

#progress bar
@app.get("/progress")
def get_progress():
    import RayTrace

    return {
        "progress": round(RayTrace.progress)
    }


@app.post("/stop-simulation")
def stop_simulation():

    import RayTrace

    RayTrace.stop_request = True

    return {
        "message": "Stop requested"
    }

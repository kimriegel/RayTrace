from fastapi import FastAPI, File, UploadFile, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, ValidationError
import sys
import json
import tempfile
from pathlib import Path
import shutil


sys.path.insert(0, "Env/")
sys.path.insert(0, "input")
sys.path.insert(0, "venv/lib/python3.7/site-packages")
sys.path.insert(0, "RayTrace/")
print(sys.path)

import newparameterfile as nPf
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


# @app.post("/run-simulation")
# def run_simulation(params: SimulationParameters):
@app.post("/run-simulation")
def run_simulation(
    params: str = Form(...),
    environmentFile: UploadFile = File(None),
    receiverFile: UploadFile = File(None),
    boomFile: UploadFile = File(None)
):
    try:
        params = SimulationParameters(**json.loads(params))
    except (ValueError, TypeError, ValidationError):
        raise HTTPException(
            status_code=422,
            detail="Invalid simulation parameters."
        )

    
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


    #NEW
    # Save the original default file paths
    default_environment = nPf.ipname
    default_receiver = nPf.RecInput
    default_boom = nPf.INPUTFILE

    # Create a temporary folder for uploaded files
    with tempfile.TemporaryDirectory() as temp_dir:

        # Start with the default files
        environment_path = default_environment
        receiver_path = default_receiver
        boom_path = default_boom

        # Save a custom environment if one was uploaded
        if environmentFile is not None:
            if Path(environmentFile.filename or "").suffix.lower() != ".obj":
                raise HTTPException(400, "Environment must be an .obj file.")

            environment_path = str(Path(temp_dir) / "environment.obj")
            with open(environment_path, "wb") as file:
                shutil.copyfileobj(environmentFile.file, file)

        # Save a custom receiver file if one was uploaded
        if receiverFile is not None:
            if Path(receiverFile.filename or "").suffix.lower() != ".txt":
                raise HTTPException(400, "Receiver must be a .txt file.")

            receiver_path = str(Path(temp_dir) / "receivers.txt")
            with open(receiver_path, "wb") as file:
                shutil.copyfileobj(receiverFile.file, file)

        # Save a custom input boom file if one was uploaded
        if boomFile is not None:
            if Path(boomFile.filename or "").suffix.lower() != ".txt":
                raise HTTPException(400, "Input boom must be a .txt file.")

            boom_path = str(Path(temp_dir) / "boom.txt")
            with open(boom_path, "wb") as file:
                shutil.copyfileobj(boomFile.file, file)

        print("Environment:", environment_path)
        print("Receiver:", receiver_path)
        print("Input boom:", boom_path)


        import RayTrace

        try:
            # Tell Python which files to use
            nPf.RecInput = receiver_path
            nPf.INPUTFILE = boom_path

            # Load the selected environment
            RayTrace.Gp.load_geometry(environment_path)

            # Run the simulation
            graph_data = RayTrace.main()

        finally:
            # Restore the original file paths
            nPf.RecInput = default_receiver
            nPf.INPUTFILE = default_boom

            # Restore the original environment
            RayTrace.Gp.load_geometry(default_environment)
    #NEW


    # import RayTrace

    # print(RayTrace)
    # print(RayTrace.__file__)
    # print(dir(RayTrace))

    # graph_data = RayTrace.main()   
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

const runBtn = document.getElementById("run-btn");
const resultParagraph = document.getElementById("result");

// progress bar add
const progressBar = document.getElementById("progress-bar");
const stopBtn = document.getElementById("stop-btn");
const progressPercent = document.getElementById("progress-percent");

stopBtn.addEventListener("click", async () => {

  await fetch("http://127.0.0.1:8000/stop-simulation", {
    method: "POST"
  });

  console.log("Stop requested");
});

runBtn.addEventListener("click", async () => {
  //progress bar
  progressBar.value = 0;
  progressPercent.textContent = "0%";

  stopBtn.disabled = false;
  runBtn.disabled = true;

  const params = {
    Fs: Number(document.getElementById("Fs").value),
    xinitial: Number(document.getElementById("xinitial").value),
    yinitial: Number(document.getElementById("yinitial").value),
    zinitial: Number(document.getElementById("zinitial").value),
    radius: Number(document.getElementById("radius").value),
    soundspeed: Number(document.getElementById("soundspeed").value),
    ps: Number(document.getElementById("ps").value),
    Temp: Number(document.getElementById("Temp").value),
    hr: Number(document.getElementById("hr").value),
    theta: Number(document.getElementById("theta").value),
    phi: Number(document.getElementById("phi").value),
    boomspacing: Number(document.getElementById("boomspacing").value),
    h: Number(document.getElementById("h").value),
  };

  try {

    // this checks python for simulation progress every half second
    const progressChecker = setInterval(async () => {

      const progressResponse = await fetch("http://127.0.0.1:8000/progress");

      const progressData = await progressResponse.json();

      console.log("Progress:", progressData.progress + "%");

      progressBar.value = progressData.progress;
      progressPercent.textContent = progressData.progress + "%";

    }, 500);


    const response = await fetch("http://127.0.0.1:8000/run-simulation", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(params),
    });



    const data = await response.json();

    clearInterval(progressChecker);

    stopBtn.disabled = true;
    runBtn.disabled = false;

    if (data.stopped) {

      // Keep the progress bar at the percentage where it stopped
      const progressResponse = await fetch("http://127.0.0.1:8000/progress");
      const progressData = await progressResponse.json();

      progressBar.value = progressData.progress;
      progressPercent.textContent = progressData.progress + "%";

    } else {

      // Only show 100% if the simulation finished normally
      progressBar.value = 100;
      progressPercent.textContent = "100%";
    }
//end


// trying to fix for stop sim button so moved rn to test 
    // console.log(data.graph_data) //new for graphs
    // console.log("Number of receivers:", data.graph_data.receivers.length); //TO CHECK IF TOO MANY FROM PYTHON

    // resultParagraph.textContent = data.message;
    resultParagraph.innerHTML = `
      <strong>${data.message}</strong><br><br>

      Fs = ${data.parameters.Fs}<br>
      xinitial = ${data.parameters.xinitial}<br>
      yinitial = ${data.parameters.yinitial}<br>
      zinitial= ${data.parameters.zinitial}<br>
      radius= ${data.parameters.radius}<br>
      soundspeed= ${data.parameters.soundspeed}<br>
      ps= ${data.parameters.ps}<br>
      Temp= ${data.parameters.Temp}<br>
      hr= ${data.parameters.hr}<br>
      theta= ${data.parameters.theta}<br>
      phi= ${data.parameters.phi}<br>
      boomspacing = ${data.parameters.boomspacing}<br>
      h = ${data.parameters.h}
    `;


    // FOR STOP SIM BUTTON
    if (data.stopped) {
      // Remove graphs from any previous simulation
      document.getElementById("graphs").innerHTML = "";

      return;
    }

    // checkiung to see if should go here 
    console.log(data.graph_data);
    console.log("Number of receivers:", data.graph_data.receivers.length);

    // gets the graph data from Python
    const graphData = data.graph_data;

    // finds HTML container for all graphs
    const graphsContainer = document.getElementById("graphs");

    // clears any graphs from a previous simulation just in case didnt refresh webpage will discard old and create new, otherwise would keep adding 5 graphs 
    graphsContainer.innerHTML = "";

    // makes a separate graph for each receiver that python sends
    graphData.receivers.forEach((receiver) => {

      // creates a new div for this receiver's graph
      const graphDiv = document.createElement("div");

      // Add the new div to the webpage
      graphsContainer.appendChild(graphDiv);


      // creates the graphs
      Plotly.newPlot(graphDiv, [
        {
          x: graphData.time,
          y: receiver.pressure,
          type: "scatter",
          mode: "lines",
          name: "Receiver " + receiver.receiver
        }
      ], {
        title: "Pressure vs Time of Receiver " + receiver.receiver,
        xaxis: {
          title: "Time [s]"
        },
        yaxis: {
          title: "Pressure [Pa]"
        }
      }, {
        responsive: true
      });

    });

  } catch (error) {
    resultParagraph.textContent = "Error connecting to Python/FastAPI.";
    console.error(error);
  }
});

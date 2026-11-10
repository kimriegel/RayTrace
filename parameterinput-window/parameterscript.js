const runBtn = document.getElementById("run-btn");
const resultParagraph = document.getElementById("result");


// GRAPH DASHBOARD ELEMENTS
const graphSection = document.getElementById("graph-section");
const receiverSelect = document.getElementById("receiver-select");
const graphsContainer = document.getElementById("graphs");

// FOR COMPARISON ELEMENTS
const singleModeBtn = document.getElementById("single-mode");
const compareModeBtn = document.getElementById("compare-mode");

const singleControls = document.getElementById("single-controls");
const compareControls = document.getElementById("compare-controls");

const receiverCheckboxes = document.getElementById("receiver-checkboxes");

// COMPARISON SELECTION BUTTONS
const selectAllBtn = document.getElementById("select-all-btn");
const resetSelectionBtn = document.getElementById("reset-selection-btn");

// Tracks which graph mode the user is viewing
let graphMode = "single";

// Stores the results from the most recent simulation
let graphData = null;

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

  // Clear and hide graphs from the previous simulation for new ones
  graphSection.style.display = "none";
  graphData = null;
  Plotly.purge(graphsContainer);


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

// ADDING FOR FILE IMPLEMENTATION
// Get the files selected by the user
    const environmentFile = document.getElementById("environmentFile").files[0];
    const receiverFile = document.getElementById("receiverFile").files[0];
    const boomFile = document.getElementById("boomFile").files[0];

    // Create a package to send parameters and files to Python
    const formData = new FormData();

    // Add all 13 simulation parameters as JSON
    formData.append("params", JSON.stringify(params));

    // Only add a file if the user selected one
    if (environmentFile) {
      formData.append("environmentFile", environmentFile);
    }

    if (receiverFile) {
      formData.append("receiverFile", receiverFile);
    }

    if (boomFile) {
      formData.append("boomFile", boomFile);
    }

    // END OF FILE IMPLEMENTATION 

let progressChecker;

  try {

    // this checks python for simulation progress every half second
    progressChecker = setInterval(async () => {

      const progressResponse = await fetch("http://127.0.0.1:8000/progress");

      const progressData = await progressResponse.json();

      console.log("Progress:", progressData.progress + "%");

      progressBar.value = progressData.progress;
      progressPercent.textContent = progressData.progress + "%";

    }, 500);


    // const response = await fetch("http://127.0.0.1:8000/run-simulation", {
    //   method: "POST",
    //   headers: {
    //     "Content-Type": "application/json",
    //   },
    //   body: JSON.stringify(params),
    // });

// THE ABOVE SECTION BEING REPLACED WITH THE BELOW FOR FILE IMPLEMENTATION
    const response = await fetch("http://127.0.0.1:8000/run-simulation", {
      method: "POST",
      body: formData,
    });

    const data = await response.json();
    
    // check if good
    if (!response.ok) {
      throw new Error(data.detail || "Simulation failed.");
    }

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

      Plotly.purge(graphsContainer);
      graphSection.style.display = "none";
      graphData = null;

      return;
    }

    // checkiung to see if should go here 
    console.log(data.graph_data);
    console.log("Number of receivers:", data.graph_data.receivers.length);

    // Store the graph data returned by Python
    graphData = data.graph_data;

    // Clear any graph from the previous simulation
    Plotly.purge(graphsContainer);

    // Clear the old dropdown options
    receiverSelect.innerHTML = "";

    // Create one dropdown option for each receiver
    graphData.receivers.forEach((receiver, index) => {

      const option = document.createElement("option");

      option.value = index;
      option.textContent = "Receiver " + receiver.receiver;

      receiverSelect.appendChild(option);

    });


    // Clear checkboxes from any previous simulation
    receiverCheckboxes.innerHTML = "";

    // Create a checkbox for every receiver
    graphData.receivers.forEach((receiver, index) => {

      const label = document.createElement("label");

      const checkbox = document.createElement("input");

      checkbox.type = "checkbox";
      checkbox.value = index;

      // Select the first two receivers automatically
      checkbox.checked = index < 2;

      // Update the graph whenever a checkbox changes
      checkbox.addEventListener("change", compareReceivers);

      label.appendChild(checkbox);

      label.appendChild(
        document.createTextNode(" Receiver " + receiver.receiver)
      );

      receiverCheckboxes.appendChild(label);

    });



    // Start each simulation in Single Receiver mode
    graphMode = "single";

    // Highlight the Single Receiver button
    singleModeBtn.classList.add("active");
    compareModeBtn.classList.remove("active");

    singleControls.style.display = "block";
    compareControls.style.display = "none";

    // Show the graph dashboard
    graphSection.style.display = "block";

    // Display the first receiver automatically
    receiverSelect.value = "0";
    showReceiver(0);

    
  } catch (error) {
    resultParagraph.textContent = "Simulation error: " + error.message;
    console.error(error);

    clearInterval(progressChecker);

    // Re-enable the Run button if something goes wrong
    runBtn.disabled = false;
    stopBtn.disabled = true;
  }

});


// DISPLAY ONE RECEIVER'S GRAPH
function showReceiver(index) {

  if (!graphData) return;

  // Get the selected receiver
  const receiver = graphData.receivers[index];

  // Create the waveform
  const trace = {
    x: graphData.time,
    y: receiver.pressure,
    type: "scatter",
    mode: "lines",
    name: "Receiver " + receiver.receiver,
    line: {
      color: "#2563eb",
      width: 2
    }
  };

  // Graph appearance
  const layout = {
    title: {
      text: "Pressure vs. Time",
      font: {
        size: 18,
        color: "#171717"
      }
    },

    xaxis: {
      title: "Time [s]",
      showgrid: true,
      gridcolor: "#e5e7eb",
      zeroline: false
    },

    yaxis: {
      title: "Pressure [Pa]",
      showgrid: true,
      gridcolor: "#e5e7eb",
      zeroline: false
    },

    hovermode: "x unified",

    showlegend: true,

    legend: {
      orientation: "h",
      x: 0.5,
      xanchor: "center",
      y: -0.3
    },

    paper_bgcolor: "white",
    plot_bgcolor: "white",

    margin: {
      l: 65,
      r: 25,
      t: 65,
      b: 100
    }
  };

  // Draw the graph PLOTLY GRAPH CONTROLS
  Plotly.react(graphsContainer, [trace], layout, {
    responsive: true,
    displaylogo: false,

    // Export graph as a PNG image
    toImageButtonOptions: {
      format: "png",
      filename: "sonic-boom-receiver-" + receiver.receiver,
      height: 650,
      width: 1100,
      scale: 2
    },

    // Remove tools we don't need
    modeBarButtonsToRemove: [
      "lasso2d",
      "select2d"
    ]
  });

}

// When the user selects another receiver
receiverSelect.addEventListener("change", () => {

  const selectedIndex = Number(receiverSelect.value);

  showReceiver(selectedIndex);

});

// THIS IS TO DISPLAY MULTIPLE RECEIVERS ON ONE GRAPH
function compareReceivers() {

  if (!graphData) return;

  // Find all the checked receiver boxes
  const checkedBoxes = document.querySelectorAll(
    '#receiver-checkboxes input[type="checkbox"]:checked'
  );

  // Create a separate trace for every selected receiver
  const traces = [];

  const colors = [
    "#2563eb", // Blue
    "#e76f51", // Orange
    "#16a085", // Teal
    "#a855f7", // Purple
    "#d69e23"  // Gold
  ];

  checkedBoxes.forEach((checkbox) => {

    // Find the receiver using its position in the list
    const index = Number(checkbox.value);

    const receiver = graphData.receivers[index];

    // Create one waveform for this receiver
    traces.push({
      x: graphData.time,
      y: receiver.pressure,
      type: "scatter",
      mode: "lines",
      name: "Receiver " + receiver.receiver,
      line: {
        color: colors[index % colors.length],
        width: 2
      }
    });

  });

  // Graph appearance
  const layout = {
    title: {
      text: "Pressure vs. Time",
      font: {
        size: 18,
        color: "#171717"
      }
    },

    xaxis: {
      title: "Time [s]",
      showgrid: true,
      gridcolor: "#e5e7eb",
      zeroline: false
    },

    yaxis: {
      title: "Pressure [Pa]",
      showgrid: true,
      gridcolor: "#e5e7eb",
      zeroline: false
    },

    hovermode: "x unified",

    showlegend: true,

    legend: {
      orientation: "h",
      x: 0.5,
      xanchor: "center",
      y: -0.3
    },

    paper_bgcolor: "white",
    plot_bgcolor: "white",

    margin: {
      l: 65,
      r: 25,
      t: 65,
      b: 100
    }
  };

  // Update the existing graph with the selected receivers
  Plotly.react(graphsContainer, traces, layout, {
    responsive: true,
    displaylogo: false,

    toImageButtonOptions: {
      format: "png",
      filename: "sonic-boom-receiver-comparison",
      height: 650,
      width: 1100,
      scale: 2
    },

    modeBarButtonsToRemove: [
      "lasso2d",
      "select2d"
    ]
  });

}


// SINGLE RECEIVER BUTTON
singleModeBtn.addEventListener("click", () => {

  graphMode = "single";

  // Show dropdown and hide checkboxes
  singleControls.style.display = "block";
  compareControls.style.display = "none";

  // Display the receiver selected in the dropdown
  showReceiver(Number(receiverSelect.value));

  singleModeBtn.classList.add("active");
  compareModeBtn.classList.remove("active");

});


// COMPARE RECEIVERS BUTTON
compareModeBtn.addEventListener("click", () => {

  graphMode = "compare";

  // Hide dropdown and show checkboxes
  singleControls.style.display = "none";
  compareControls.style.display = "block";

  // Display all checked receivers together
  compareReceivers();

  compareModeBtn.classList.add("active");
  singleModeBtn.classList.remove("active");

});


// SELECT ALL RECEIVERS
selectAllBtn.addEventListener("click", () => {

  const checkboxes = receiverCheckboxes.querySelectorAll(
    'input[type="checkbox"]'
  );

  // Check every receiver
  checkboxes.forEach((checkbox) => {
    checkbox.checked = true;
  });

  // Update the graph
  compareReceivers();

});


// RESET RECEIVER SELECTION
resetSelectionBtn.addEventListener("click", () => {

  const checkboxes = receiverCheckboxes.querySelectorAll(
    'input[type="checkbox"]'
  );

  // Select only the first two receivers
  checkboxes.forEach((checkbox, index) => {
    checkbox.checked = index < 2;
  });

  // Update the graph
  compareReceivers();

});
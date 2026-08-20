const runBtn = document.getElementById("run-btn");
const resultParagraph = document.getElementById("result");

runBtn.addEventListener("click", async () => {
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

    const response = await fetch("http://127.0.0.1:8000/run-simulation", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(params),
    });

    const data = await response.json();

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
  } catch (error) {
    resultParagraph.textContent = "Error connecting to Python/FastAPI.";
    console.error(error);
  }
});


// const sendBtn = document.getElementById("send-btn");

// sendBtn.addEventListener("click", async () => {
//   const userInput = document.getElementById("user-input").value;

//   const response = await fetch(
//     `http://127.0.0.1:8000/inputparameters?h=${encodeURIComponent(userInput)}`
//   );

//   const data = await response.json();

//   resultParagraph.textContent = data.message;
// });
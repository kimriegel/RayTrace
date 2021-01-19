//Interact with trace file

// Update bar
const progressBar = document.getElementsByClassName('progress-bar')[0]
setInterval(() => {
    const computedStyle = getComputedStyle(progressBar)
    const width = parseFloat(computedStyle.getPropertyValue('--width')) || 0
    progressBar.style.setProperty('--width',width + .1)
    // const width = rayNo/RAYMAX *100
    //progressBar.style.setProperty('--width',width)
    
}, 5)
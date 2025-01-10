// Purpose: JavaScript code for the interactive features.
document.addEventListener("DOMContentLoaded", () => {
    const speechToSmilesButton = document.getElementById("speechToSmiles"); // Button to start voice recognition
    const smilesInput = document.getElementById("smiles"); // Input field for SMILES strings

    // Drug names to their corresponding SMILES strings
    const drugToSmilesMap = {
        aspirin: "CC(=O)OC1=CC=CC=C1C(=O)O",
        paracetamol: "CC(=O)NC1=CC=C(C=C1)O",
        ibuprofen: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        glucose: "C6H12O6",
    };

    // Event listener for the voice recognition button
    speechToSmilesButton.addEventListener("click", () => {
        const SpeechRecognition =
            window.SpeechRecognition || window.webkitSpeechRecognition; // Cross-browser support for SpeechRecognition
        const recognition = new SpeechRecognition();

        // Event triggered when speech is successfully recognized
        recognition.onresult = (event) => {
            const transcript = event.results[0][0].transcript.trim().toLowerCase(); // Capture the recognized text
            console.log("You said:", transcript);

            // If the recognized text matches a known drug name, convert it to SMILES
            if (drugToSmilesMap[transcript]) {
                smilesInput.value = drugToSmilesMap[transcript]; // Populate the input with the SMILES string
                console.log(`Converted ${transcript} to SMILES: ${drugToSmilesMap[transcript]}`);
            } else {
                // If not a known drug, assume it's a SMILES string
                smilesInput.value = transcript;
                console.log("Input recognized as SMILES:", transcript);
            }
        };

        // Handle any errors during speech recognition
        recognition.onerror = (event) => {
            console.error("Speech recognition error:", event.error);
            alert("Speech recognition failed. Please try again.");
        };

        recognition.start(); // Start listening for speech
    });
});

// Drawing tool functionality
document.addEventListener("DOMContentLoaded", () => {
    const canvas = document.getElementById("molecule-drawing"); // Canvas element for drawing
    const ctx = canvas.getContext("2d"); // Get the 2D drawing context
    let isDrawing = false; // Track whether the user is currently drawing

    // Set drawing properties
    ctx.lineWidth = 2; // Set the thickness of the drawing line

    // Start drawing when the mouse is pressed down
    canvas.addEventListener("mousedown", (e) => {
        isDrawing = true;
        ctx.beginPath(); // Start a new path
        ctx.moveTo(e.offsetX, e.offsetY); // Move to the position of the mouse
    });

    // Draw as the mouse moves
    canvas.addEventListener("mousemove", (e) => {
        if (isDrawing) {
            ctx.lineTo(e.offsetX, e.offsetY); // Draw a line to the current mouse position
            ctx.stroke(); // Render the line on the canvas
        }
    });

    // Stop drawing when the mouse is released
    canvas.addEventListener("mouseup", () => {
        isDrawing = false;
    });

    // Stop drawing if the mouse leaves the canvas
    canvas.addEventListener("mouseout", () => {
        isDrawing = false;
    });

    // Clear the canvas when the clear button is clicked
    document.getElementById("clearCanvas").addEventListener("click", () => {
        ctx.clearRect(0, 0, canvas.width, canvas.height); // Clears the canvas
    });
});

// Placeholder for search button functionality
document.getElementById('search-button').addEventListener('click', function() {
    alert('Search functionality is not implemented yet.'); // Pending message for search functionality
});


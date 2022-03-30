<!DOCTYPE html>
<html lang="en-GB">
<head>
    <meta charset="UTF-8">
    <title>Graph</title>
    <link rel="stylesheet" href="media/styles.css">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Atkinson+Hyperlegible:ital,wght@0,400;0,700;1,400;1,700&display=swap" rel="stylesheet">
</head>
<body class="graph">
<?php

// function to redirect the page to a different url
function redirect($url) {
    ob_start();
    header('Location: '.$url);
    ob_end_flush();
    die();
}

// extract parameters from the form
$latex     = $_GET["latex"];
$lowerx    = $_GET["lowerx"];
$upperx    = $_GET["upperx"];
$lowery    = $_GET["lowery"];
$uppery    = $_GET["uppery"];
$linestep  = $_GET["linestep"];
$precision = $_GET["precision"];
$centerx   = $_GET["centerx"];
$centery   = $_GET["centery"];
$radius    = $_GET["radius"];
$clientw   = $_GET["clientw"];
$clienth   = $_GET["clienth"];

// set time limit on execution of 2 hours
set_time_limit(7200);

// create command with appropriate parameters
$command = "autosig\scripts\python.exe main.py \"{$latex}\" {$lowerx} {$upperx} {$lowery} {$uppery} {$linestep} {$precision} {$centerx} {$centery} {$radius} {$clientw} {$clienth}";

echo $command;

// execute command
$output = shell_exec($command);

// echo $output;

// get the raw query string from the server
$query = $_SERVER['QUERY_STRING'];

// if no output was given, that means an unknown error has been encountered
if (!$output) {
    $output = "EUnknown Error";
}

// if an error has been encountered, then send the user back to index.php with the error message as an argument
if ($output[0] == "E") {
    $error = rawurlencode(substr($output, 1));
    redirect("index.php?error={$error}&{$query}");
}

?>
<canvas id="lineCanvas" class="canvas" width="<?php echo $clientw;?>" height="<?php echo $clienth;?>"></canvas>
<canvas id="gridCanvas" class="canvas" width="<?php echo $clientw;?>" height="<?php echo $clienth;?>"></canvas>
<div class="slideBox">
    <input type="range" min="0" max="1000" value="0" class="slider" id="time" oninput="update()" disabled="true"> <br>
    <button id="back" class="button" onclick="goBack()">Return To input</button>
    <button id="replay" class="button" onclick="replay()">Replay Animation</button>
</div>
<script>

// get canvas contexts and objects
lineCanvas = document.getElementById("lineCanvas")
lCtx = lineCanvas.getContext("2d");
gridCanvas = document.getElementById("gridCanvas")
gCtx = gridCanvas.getContext("2d");

// set the line canvas's linewidth to 2
lCtx.lineWidth = 2;

// set the grid canvas's alpha to 0.3
gCtx.globalAlpha = 0.3

// this works for some reason
lines = <?php echo $output;?>;

// function which produces an ease-in-out curve
function ease(t) {
    t /= 1000;
    if (t > 0.5) {
        return 4 * Math.pow(t - 1, 3) + 1;
    } else {
        return 4 * Math.pow(t, 3);
    }
}

// given a point (which actually contains an input point and an output point) and a time, return a point directly
// between those 2 points, with a position corresponding to the ease-in-out curve with respect to time
function interpolate(point, t) {
    let factor = ease(t);
    let antiFactor = 1 - factor;
    let x = factor * point[2] + antiFactor * point[0];
    let y = factor * point[3] + antiFactor * point[1];
    return [x, y]
}

// draw a singular line part. these are contiguous
function draw_linepart(points, t, context) {
    // calculate the first point
    let first = interpolate(points[0], t);

    // initialise the path
    context.beginPath();
    context.moveTo(first[0], first[1]);

    // since the points list is interweaved, then this loop must be able to unweave it
    for (let i = 0; i < (points.length - 1) / 2; i++) {
        let cntrl = points[2 * i + 1];
        let next  = points[2 * i + 2];

        // find the appropriate point inbetween the initial point and the final point
        let cpoint = interpolate(cntrl, t);
        let npoint = interpolate(next, t);

        // draw a quadratic bezier curve from the last endpoint to the next endpoint
        context.quadraticCurveTo(cpoint[0], cpoint[1], npoint[0], npoint[1]);
    }

    // actually draw the line
    context.stroke()
}

// draws all the lines at a specific point in time
function draw(lines, t, context) {

    // clear the canvas
    context.clearRect(0, 0, lineCanvas.width, lineCanvas.height);

    // draw horizontal lines
    let hor = lines["horizontal"];

    // set the line colour to a colourblind friendly red
    // source for this colour: https://davidmathlogic.com/colorblind/#%23005AB5-%23DC3220
    context.strokeStyle = '#DC3220';

    // loop through each line
    for (let i = 0; i < hor.length; i++) {
        let line = hor[i];

        // loop through each line part in this line
        for (let j = 0; j < line.length; j++) {
            draw_linepart(line[j], t, context)
        }
    }

    // draw vertical lines
    let vert = lines["vertical"];

    // set the line colour to a colourblind friendly blue
    // source for this colour: https://davidmathlogic.com/colorblind/#%23005AB5-%23DC3220
    context.strokeStyle = '#005AB5';

    // loop though each line
    for (let i = 0; i < vert.length; i++) {
        let line = vert[i];

        // loop though each line part in this line
        for (let j = 0; j < line.length; j++) {
            draw_linepart(line[j], t, context)
        }
    }
}

// when page is loaded, start drawing and animating the lines
window.onload = function() {

    // draw inital gridlines
    draw(lines, 0, gCtx);

    // draw outputted lines before any transformation
    draw(lines, 0, lCtx);
    replay();
}

// function to play/replay the animation
function replay() {

    // disable the replay button and time slider
    document.getElementById("time").disabled = true;
    document.getElementById("replay").disabled = true;
    let slider = document.getElementById("time");

    // loop though every timestep
    for (let t = 0; t <= 1000; t++) {

        // setTimeout is an asynchronous function, meaning it doesnt execute immediately, but instead allows the code
        // after it to execute while it is still executing
        // this means, to get these operations to happen in order, the timeout must be increased by the timestep
        // 500ms is also added so that it doesnt animate immediately, but gives the user time to focus on the window
        setTimeout(() => {
            draw(lines, t, lCtx);
            slider.value = t;
        }, 10*t + 500);
    }

    // do final step, which also includes re-enabling the replay button and time slider
    setTimeout(() => {
            draw(lines, 1000, lCtx); slider.value = 1000;
            document.getElementById("time").disabled = false;
            document.getElementById("replay").disabled = false;
        }, 10500);

}

// update drawing to the current value of the slider
function update() {
    let time = document.getElementById("time").value;
    console.log(time);
    draw(lines, time, lCtx);
}

// allow the user to return to the index page with the same parameters
function goBack() {
    let query = "<?php echo $_SERVER['QUERY_STRING'];?>";
    window.location.href = "index.php?" + query;
}

</script>

</body>
</html>
var canvas = document.getElementById("GameCanvas");
var ctx = canvas.getContext("2d");

/*
coord: {x: ..., y: ...}
e.g. testcoord = {x: 4, y: 6}
testcoord.x = 12
circpos: {start: ..., end: ...}
*/

// resize for mobile devices

canvas.width  = window.innerWidth;
canvas.height = window.innerHeight;

function debug_draw_text(str) {
    ctx.font = "30px Arial";
    ctx.fillText(str, 10, 50);
}

function refresh_canvas(color) {
    ctx.beginPath();
    ctx.fillStyle = color;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.closePath();
}

function draw_line(coords, thickness) {
    ctx.beginPath();
    ctx.moveTo(coords[0].x, coords[0].y);
    ctx.lineTo(coords[1].x, coords[1].y);
    ctx.lineWidth = thickness;
    ctx.stroke(); 
}

function draw_circle_filled(coord, radius, color) {
    ctx.beginPath();
    ctx.arc(coord.x, coord.y, radius, 0, 2 * Math.PI);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.closePath();
}

function draw_circle_segment(coord, radius, circpos, thickness) {
    ctx.beginPath();
    ctx.arc(coord.x, coord.y, radius, circpos.start, circpos.end);
    ctx.lineWidth = thickness;
    ctx.stroke();
    ctx.closePath();
}
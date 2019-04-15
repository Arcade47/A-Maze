// game principle
/*
circular maze, a central circle (player) needs to get out in time
problem: circle grows - when too large, player dies (crushed between maze walls)
circle size can be resetted when collecting items ("meds")
player control via arrow keys
---
conceptual TODO:
- how can mazes be generated automatically?
*/

// functions
function get_exact_coord(angle, dist, in_deg) {
    // trigonometric functions in JS are radians-based
    function deg_to_rad(deg) {
        return deg*(Math.PI/180);
    }
    if (in_deg) {
        var rad = deg_to_rad(angle); // convert
    } else {
        var rad = angle; // no change necessary
    }
    var x_val = dist*Math.cos(rad) + canvas.width/2;
    var y_val = dist*Math.sin(rad) + canvas.height/2;
    return {x: x_val, y: y_val};
}

function random_wall_rads(n_walls) {
    // returns position along ring in radians for n walls (equally spaced)
    var wall_rads = [];
    var step = (2*Math.PI)/n_walls;
    var start_rad = Math.random()*(2*Math.PI);
    for (let index = 0; index < n_walls; index++) {
        wall_rads.push(start_rad);
        start_rad += step;
        start_rad = start_rad%(2*Math.PI);
    }
    return wall_rads;
}

function overlap_hole_region(wall_rad, coords) {

    // TODO test whether wrap around is correct

    for (let j = 0; j < coords.length-1; j++) {
        var hole_start = coords[j].end;
        var hole_end = coords[j+1].start;
        // consider case that whole is around radians 0
        if (hole_end < hole_start) {
            hole_end += Math.PI*2;
            wall_rad += Math.PI*2;
        }
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            return true;
        }
    }
    // test last hole
    if (coords.length > 0) {
        var hole_start = coords[coords.length-1].end;
        var hole_end = coords[0].start;
        // consider case that whole is around radians 0
        if (hole_end < hole_start) {
            hole_end += Math.PI*2;
            wall_rad %= Math.PI*2; // remove the first addition if applied
            wall_rad += Math.PI*2;
        }
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            return true;
        }
    }

    // no overlap
    return false;

}

function walls_in_holes(coords1, coords2, wall_rads) {
    for (let i = 0; i < wall_rads.length; i++) {
        const wall_rad = wall_rads[i];
        // if any overlap with hole of either inner or outer ring --> function returns true
        if (overlap_hole_region(wall_rad, coords1)) { return true; }
        if (overlap_hole_region(wall_rad, coords2)) { return true; }
    }
    // no overlap
    return false;
}

function get_wall_validity(coords1, coords2, wall_rads) {
    var wall_list = [];
    for (let i = 0; i < wall_rads.length; i++) {
        const wall_rad = wall_rads[i];
        if (overlap_hole_region(wall_rad, coords1)) { 
            wall_list.push(false);
        }
        else if (overlap_hole_region(wall_rad, coords2)) {
            wall_list.push(false);
        }
        else {
            wall_list.push(true);
        }
    }
    return wall_list;
}

function circumference_coords(n_segs, r, ring_space, thickness, rot) {
    // get distance of start and end pos based on radius and ring_space
    // 1. calculate circumference
    var c = 2*Math.PI*r;
    // 2. standardize the ring_space
    var standard_ring_space = (ring_space + thickness)/c; // leave room for circular edges          thickness
    // 3. add correct measure for canvas arc method (2*pi)
    var end_start_dist = standard_ring_space*2*Math.PI;

    // get coordinates based on number of segments
    var step = 2*Math.PI/n_segs - end_start_dist; // - (n_segs*end_start_dist);
    // apply rotation (in percent of max)
    var rot_add = rot*2*Math.PI;
    var current_start = (end_start_dist/2 + rot_add)%(2*Math.PI); // so that first hole is centered
    var current_end = (current_start + step)%(2*Math.PI);
    var output_edges = [];          // for drawing the circular segments
    var output_ball_coords = [];    // for drawing the circular edges of the segments
    for (let index = 0; index < n_segs; index++) {
        // append to segment list
        output_edges.push({start: current_start, end: current_end});
        // append to edge coord list
        output_ball_coords.push(get_exact_coord(current_start, r, false));
        output_ball_coords.push(get_exact_coord(current_end, r, false));
        // step coords forward
        current_start = (current_end + end_start_dist)%(2*Math.PI);
        current_end = (current_start + step)%(2*Math.PI);
    }
    return [output_edges, output_ball_coords];
}

// init global vars
var n_rings;
var ring_spacing;
var ring_thickness;
var n_meds;
var size_grow;
var player_size;
var player_speed;
var max_distance;
var med_min_distance;
var med_max_distance;
var rings = [];
var walls = [];

// add event listeners
document.addEventListener("keydown", keydown);
document.addEventListener("keyup", keyup);

// prepare values for vars
function add_values() {
    // global vars that are initialized before
    n_rings = 10;            // rings of maze
    ring_spacing = 20;      // spacing of maze rings in px
    ring_thickness = 4;    // diameter of maze walls in px
    n_meds = 5;             // how many items need to be collected before allowed to exit
    size_grow = 0.005;          // increase of player diameter over a frame
    player_size = 1;        // initial player diameter in px
    player_speed = 3;       // movement speed 
    med_min_distance_fac = 0.2;  // factor of max distance that player is separated by med
    med_max_distance_fac = 0.8;  // how far is last med from player?
    // derived vars
    // TODO: max distance   // formula that takes into account player_speed, size_grow and ring_spacing
    // max_distance = ;
}

// classes

class Player {
    constructor() {
        // appearance
        this.pos = {x: canvas.width/2, y: canvas.height/2};
        this.speed = player_speed;
        this.radius = player_size;
        this.color = "black";
        // states
        this.leftmove = false;
        this.rightmove = false;
        this.upmove = false;
        this.downmove = false;
    }
    update() {
        // apply movements
        if (this.leftmove) { this.pos.x -= player_speed; }
        if (this.rightmove) { this.pos.x += player_speed; }
        if (this.upmove) { this.pos.y -= player_speed; }
        if (this.downmove) { this.pos.y += player_speed; }
        // increase in size
        // this.radius += size_grow;
    }
    render() {
        draw_circ(this.radius, this.pos, this.color);
    }
}

class Wall {
    constructor(ind1, ind2, rad, invalid=false) {
        this.ind1 = ind1; // inner ring to which wall is connected
        this.ind2 = ind2; // outer ring to which wall is connected
        this.rad = rad; // angle in radians
        // derived vars
        this.r1 = ((ind1+1)*2*ring_spacing - ring_thickness)/2;
        this.r2 = ((ind2+1)*2*ring_spacing - ring_thickness)/2;
        this.coord1 = get_exact_coord(this.rad, this.r1, false);
        this.coord2 = get_exact_coord(this.rad, this.r2, false);
        this.invalid = invalid; // for debugging
        this.moving = false;
    }
    set_validity() {
        // get them coords
        coords1 = rings[this.ind1].holes_coords;
        coords2 = rings[this.ind2].holes_coords;
        // if any overlap with hole of either inner or outer ring --> function sets validity
        this.invalid = false;
        // console.log(coords1.length);
        if (overlap_hole_region(this.rad, coords1)) { this.invalid = true; }
        if (overlap_hole_region(this.rad, coords2)) { this.invalid = true; }
    }
    update_coords() {
        this.coord1 = get_exact_coord(this.rad, this.r1, false);
        this.coord2 = get_exact_coord(this.rad, this.r2, false);
    }
    update() {
        this.set_validity();
        // position updates (debugging)
        if (this.moving) {
            this.rad += 0.05
            this.rad %= Math.PI*2; // make sure value stays within bounds
            this.update_coords();
        }
    }
    render() {
        var color = "black"
        if (this.invalid) {
            color = "red";
        }
        draw_line([this.coord1, this.coord2], color, ring_thickness);
    }
}

class Ring {
    constructor(ind, rotation) {
        // lower inds --> smaller rings
        this.diameter = (ind+1)*2*ring_spacing - ring_thickness;
        var coords_container = circumference_coords(ind+1, this.diameter/2, ring_spacing, ring_thickness, rotation);
        this.holes_coords = coords_container[0];
        this.edge_ball_coords = coords_container[1];
    }
    render() {
        // draw line segments
        for (let index = 0; index < this.holes_coords.length; index++) {
            draw_circ_segment(this.diameter/2, {x: canvas.width/2, y: canvas.height/2}, "black", ring_thickness, this.holes_coords[index]);
        }
        // draw edges
        for (let index = 0; index < this.edge_ball_coords.length; index++) {
            draw_circ(ring_thickness/2, this.edge_ball_coords[index], "black");
        }
    }
}

// instantiate objects
add_values();
var player = new Player();
// add rings and walls
for (let index = 0; index < n_rings; index++) {
    rings.push(new Ring(index, Math.random()));

    // add walls (TODO: better algorithm)
    if (index == 1) { // (index > 0)
        // make a test wall that moves upon enter presses
        // walls.push(new Wall(index, index-1, wall_rads[index2], validities[index2]));
        // /*
        // repeat finding random wall radians position until no overlap with holes
        var coords1 = rings[index].holes_coords;
        var coords2 = rings[index-1].holes_coords;
        // var wall_rads = random_wall_rads(index+1);
        var wall_rads = random_wall_rads(1);
        var validities = get_wall_validity(coords1, coords2, wall_rads);
        // find new random wall rads as long as wall in hole
        // while (walls_in_holes(coords1, coords2, wall_rads)) {
        //     console.log('new coords necessary');
        //     var wall_rads = random_wall_rads(index+1);
        // }
        for (let index2 = 0; index2 < wall_rads.length; index2++) {
            walls.push(new Wall(index, index-1, wall_rads[index2], validities[index2]));
        }
        // */
    }
}

// overall update function
function update() {
    // run changes (update objects)
    player.update();
    for (let index = 0; index < walls.length; index++) {
        walls[index].update();
    }
    // draw all changes
    draw();
    // get animation going
    requestAnimationFrame(update);
}

// overall draw function
function draw() {
    // refresh
    set_canvas_bg("lightblue");
    // draw maze
    for (let index = 0; index < rings.length; index++) {
        rings[index].render();
    }
    for (let index = 0; index < walls.length; index++) {
        walls[index].render();
        
    }
    // draw player
    player.render();
}

// event listener actions
function keydown(e) {
    // TODO make gameplay more fun - diagonal movements?
    // left
    if (e.keyCode == 37) { player.leftmove = true; }
    // up
    if (e.keyCode == 38) { player.upmove = true; }
    // right
    if (e.keyCode == 39) { player.rightmove = true; }
    // down
    if (e.keyCode == 40) { player.downmove = true; }

    // enter
    if (e.keyCode == 13) { 
        for (let index = 0; index < walls.length; index++) {
            walls[index].moving = true;
        }
    }
}
function keyup(e) {
    // left
    if (e.keyCode == 37) { player.leftmove = false; }
    // up
    if (e.keyCode == 38) { player.upmove = false; }
    // right
    if (e.keyCode == 39) { player.rightmove = false; }
    // down
    if (e.keyCode == 40) { player.downmove = false; }

    // enter
    if (e.keyCode == 13) { 
        for (let index = 0; index < walls.length; index++) {
            walls[index].moving = false;
        }
    }
}

// start game loop

update();
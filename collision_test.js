// game principle
/*
circular maze, a central circle (player) needs to get out in time
problem: circle grows - when too large, player dies (crushed between maze walls)
circle size can be resetted when collecting items ("meds")
player control via arrow keys
*/

// init global vars
var n_rings;
var ring_spacing;
var ring_thickness;
var n_meds;
var size_grow;
var player_size;
var player_speed;
var collision_steps;
var max_holes;
var dead_end_frac;
var max_distance;
var med_min_distance;
var med_max_distance;
var rings = [];
var walls = [];
var start_rot;
var hole_problem_level = 1; // layer at which overlap problems need to be resolved

// debug flags
var debug_draw_toggle = false;

// add event listenersggg
document.addEventListener("keydown", keydown);
document.addEventListener("keyup", keyup);
document.addEventListener("click", mouseclick);

// prepare values for vars
function add_values() {
    // global vars that are initialized before
    n_rings = 3;            // rings of maze
    ring_spacing = 80;      // spacing of maze rings in px
    ring_thickness = 20;    // diameter of maze walls in px
    n_meds = 5;             // how many items need to be collected before allowed to exit
    size_grow = 0.005;          // increase of player diameter over a frame
    player_size = 15;        // initial player diameter in px
    player_speed = 5;       // movement speed
    collision_steps = 10;   // determines how many times per draw the collision resolution function is called
    max_holes = 6;          // number of holes for middle rings
    dead_end_frac = 0.2;       // percentage of nodes in typical layer that are dead ends
    med_min_distance_fac = 0.2;  // factor of max distance that player is separated by med
    med_max_distance_fac = 0.8;  // how far is last med from player?
    start_rot = Math.random(); // rotation of innermost ring
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
        this.color = "green";
        // states
        this.leftmove = false;
        this.rightmove = false;
        this.upmove = false;
        this.downmove = false;
        // for collision
        this.vel = {x: 0, y: 0};
        this.current_ring_ind = 0; // starts in middle
        this.collision_balls = [];
    }
    set_relevant_coll_balls() {
    }
    resolve_collisions(old_pos) {
    }
    update() {
        // copy old pos (before changes) for velocity derivation
        var old_pos = {x: this.pos.x, y: this.pos.y};

        // apply movements
        if (this.leftmove) { this.pos.x -= player_speed; }
        if (this.rightmove) { this.pos.x += player_speed; }
        if (this.upmove) { this.pos.y -= player_speed; }
        if (this.downmove) { this.pos.y += player_speed; }

        // derive vel
        this.vel.x = this.pos.x - old_pos.x;
        this.vel.y = this.pos.y - old_pos.y;

        // collision resolution (in multiple steps)
        this.resolve_collisions(old_pos);

        // update virtual collision objects
        this.set_relevant_coll_balls();

        // increase in size
        // this.radius += size_grow;
    }
    render() {

        // debug: render collision_balls
        for (let index = 0; index < this.collision_balls.length; index++) {
            const cball = this.collision_balls[index];
            draw_circ(ring_thickness/2, cball, "white");
        }

        draw_circ(this.radius, this.pos, this.color);

        // debug: draw vel
        draw_line([this.pos, {x: this.pos.x + this.vel.x, y: this.pos.y + this.vel.y}], "white", 3);

    }
}

class TestGrid {
    constructor() {

    }
    render() {
        // one ring
        draw_circ_outline_transp(ring_spacing, center_coord, "black", ring_thickness);
    }
}

// instantiate objects
add_values();
var player = new Player();
var testgrid = new TestGrid();

// overall update function
function update() {
    // run changes (update objects)
    player.update();
    // draw all changes
    draw();
    // get animation going
    requestAnimationFrame(update);
}

// overall draw function
function draw() {
    // refresh
    set_canvas_bg("lightblue");
    if (debug_draw_toggle) {
    }
    // draw testgrid
    testgrid.render();
    // draw player
    player.render();
}

// event listener actions
function keydown(e) {
    // left
    if (e.keyCode == 37) { player.leftmove = true; }
    // up
    if (e.keyCode == 38) { player.upmove = true; }
    // right
    if (e.keyCode == 39) { player.rightmove = true; }
    // down
    if (e.keyCode == 40) { player.downmove = true; }

    // enter --> debug rotate walls
    if (e.keyCode == 13) {
    }

    // S, W --> debug increase/decrease player size
    if (e.keyCode == 87) { player.radius += 20*size_grow; }
    if (e.keyCode == 83 && player.radius > 1 + size_grow) { player.radius -= 20*size_grow; }
    // B --> toggle debug
    if (e.keyCode == 66) { debug_draw_toggle = !debug_draw_toggle; }

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
    }
}
function mouseclick(e) {
    // debug: select cell
    var pos = getXY(e);
    player.pos.x = pos.x;
    player.pos.y = pos.y;
}

// start game loop
update();
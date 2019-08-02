// game principle
/*
circular maze, a central circle (player) needs to get out in time
problem: circle grows - when too large, player dies (crushed between maze walls)
circle size can be resetted when collecting items ("meds")
player control via arrow keys
*/

// init global vars
var min_rings = 5;
var max_rings = 12;
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
var player;
var med; // TODO adjust that it can be multiple meds
var network;
var maze;

// add event listenersggg
document.addEventListener("keydown", keydown);
document.addEventListener("keyup", keyup);

// prepare values for vars
function reset(nrings) {
    // global vars that are initialized before
    n_rings = nrings;            // rings of maze
    ring_spacing = 30;      // spacing of maze rings in px
    ring_thickness = 5;    // diameter of maze walls in px
    n_meds = 1;//Math.floor(n_rings/2);             // how many items need to be collected before allowed to exit
    size_grow = 0.0075;          // increase of player diameter over a frame
    player_size = 1.5;        // initial player diameter in px
    player_speed = 2;       // movement speed
    collision_steps = 30;   // determines how many times per draw the collision resolution function is called
    max_holes = 6;          // number of holes for middle rings
    dead_end_frac = 0.2;       // percentage of nodes in typical layer that are dead ends
    med_min_distance_fac = 0.2;  // factor of max distance that player is separated by med
    med_max_distance_fac = 0.8;  // how far is last med from player?
    start_rot = Math.random(); // rotation of innermost ring
    // derived vars
    // TODO: max distance   // formula that takes into account player_speed, size_grow and ring_spacing
    // max_distance = ;
    
    // init classes
    player = new Player();
    network = new MazeNetworkRandom();
    maze = new MazeRandom(network);
    med = new Med();
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
        // debug movements
        this.CWmove = false;
        this.CCWmove = false;
        // for collision
        this.vel = {x: 0, y: 0};
        this.current_ring_ind = 0; // starts in middle
        this.collision_balls = [];
    }
    set_current_ring_ind() {
        // function of distance to center, pos and ring_spacing (ring_thickness?)
        var dist = distance(this.pos, center_coord);
        // set value of current ring ind
        this.current_ring_ind = Math.floor(dist/ring_spacing);
    }
    set_relevant_coll_balls_ring(rad, ring_ind, dist) {
        var relevant_holes_coords = maze.rings[ring_ind].holes_coords;
        // TODO write function that it makes sense
        var hole_overlap = overlap_hole_region(rad, relevant_holes_coords, false, false);
        if (hole_overlap) {
            var hole_coords = overlap_hole_region(rad, relevant_holes_coords, true, false)[1];
            // console.log(hole_coords);
            var coord1 = get_exact_coord(hole_coords.start, dist, false);
            var coord2 = get_exact_coord(hole_coords.end, dist, false);
            this.collision_balls.push(coord1);
            this.collision_balls.push(coord2);
            // console.log("OVERLAP");
        } else {
            var coll_ball_pos = get_exact_coord(rad, dist, false);
            this.collision_balls.push(coll_ball_pos);
            // console.log("NO OVERLAP");
        }
    }
    set_relevant_coll_balls_wall(player_dist_center) {
        for (let index = 0; index < maze.rings[this.current_ring_ind].walls.length; index++) {
            const wall = maze.rings[this.current_ring_ind].walls[index];
            // convert to coords
            if (!wall.dug) {
                this.collision_balls.push(get_exact_coord(wall.rad, player_dist_center, false));
            }
        }
    }
    set_relevant_coll_balls() {
        // deriving current ring ind
        this.set_current_ring_ind();
        // reset coll balls
        this.collision_balls = [];
        // first set coll balls along ring
        var inner_dist = this.current_ring_ind*ring_spacing - 0.5*ring_thickness;
        var outer_dist = (this.current_ring_ind + 1)*ring_spacing - 0.5*ring_thickness;
        var player_dist_center = distance(this.pos, center_coord);
        var player_rad = coord_to_rad(this.pos);
        player_rad = flip_rad(player_rad); // TODO find out why necessary
        // ignore the innermost pseudoring (-1)
        if (this.current_ring_ind > 0 && this.current_ring_ind < n_rings) {
            this.set_relevant_coll_balls_ring(player_rad, this.current_ring_ind - 1, inner_dist);
            this.set_relevant_coll_balls_wall(player_dist_center);
        }
        if (this.current_ring_ind < n_rings) {
            this.set_relevant_coll_balls_ring(player_rad, this.current_ring_ind, outer_dist);
            this.set_relevant_coll_balls_wall(player_dist_center);
        }
    }
    resolve_collisions() {

        // run in simulation steps: partition the velocity
        var vel_step = {x: this.vel.x/collision_steps, y: this.vel.y/collision_steps};
        // run the collision check collision_steps times --> only as long as no obstacle is reached
        var obstacle = false;

        // cycle through every velocity component
        for (let i = 0; i < collision_steps; i++) {
            // only update pos if no obstacle
            if (!obstacle) {
                // step forward the position
                this.pos.x += vel_step.x;
                this.pos.y += vel_step.y;

                // set the coll balls to new position
                this.set_relevant_coll_balls();

                // loop through all collision balls, check collision
                // circles thru the collision balls
                for (let j = 0; j < this.collision_balls.length; j++) {
                    const cball = this.collision_balls[j];
                    // check if overlap --> collision needs to be resolved
                    if (distance(this.pos, cball) < (this.radius + ring_thickness/2)) {
                        // get the overlap
                        var overlap = (this.radius + ring_thickness/2) - distance(this.pos, cball);
                        // get the overlap vector
                        var xo = cball.x - this.pos.x;
                        var yo = cball.y - this.pos.y;
                        var overlap_vec = {x: xo, y: yo};
                        // get the length of overlap vector
                        var l = Math.sqrt(overlap_vec.x*overlap_vec.x + overlap_vec.y*overlap_vec.y);
                        // console.log(l);
                        this.pos.x -= overlap*overlap_vec.x/l;
                        this.pos.y -= overlap*overlap_vec.y/l;
                    }
                }
            }
        }
    }
    win() {
        // player distance to center large enough?
        if (this.current_ring_ind > n_rings) {
            if (n_rings < max_rings) {
                var new_n_rings = n_rings + 1;
            } else {
                var new_n_rings = max_rings;
            }
            reset(new_n_rings);
        }
    }
    lose() {
        // cycle through all collision balls again: if remaining overlap - death
        this.set_relevant_coll_balls();
        for (let index = 0; index < this.collision_balls.length; index++) {
            const cball = this.collision_balls[index];
            if (distance(this.pos, cball) < (this.radius + ring_thickness/2)) {
                var overlap = (this.radius + ring_thickness/2) - distance(this.pos, cball);
                if (overlap > 0.1) {
                    if (n_rings == min_rings) {
                        // TODO: game over
                        reset(n_rings);
                    } else {
                        reset(n_rings - 1);
                    }
                    break;
                }
            }
        }
    }
    update() {
        // copy old pos (before changes) for velocity derivation
        var old_pos = {x: this.pos.x, y: this.pos.y};

        // apply movements
        if (this.leftmove) { this.pos.x -= player_speed; }
        if (this.rightmove) { this.pos.x += player_speed; }
        if (this.upmove) { this.pos.y -= player_speed; }
        if (this.downmove) { this.pos.y += player_speed; }

        // debug: circular movements (keys Q, E)
        if (this.CWmove) { this.circular_move(false); }
        if (this.CCWmove) { this.circular_move(true); }

        // derive vel
        this.vel.x = this.pos.x - old_pos.x;
        this.vel.y = this.pos.y - old_pos.y;

        // collision resolution (in multiple steps)
        this.resolve_collisions();

        // update current position in terms of ring indices
        this.set_current_ring_ind();

        // update virtual collision objects
        // this.set_closest_collision_balls();
        this.set_relevant_coll_balls();

        // increase in size
        this.radius += size_grow;

        // check if maze solved or player crushed
        this.win();
        this.lose();
    }
    render() {

        // debug: render collision_balls
        /*
        for (let index = 0; index < this.collision_balls.length; index++) {
            const cball = this.collision_balls[index];
            draw_circ(ring_thickness/2, cball, "white");
        }
        */

        draw_circ(this.radius, this.pos, this.color);

        // debug: draw vel
        // draw_line([this.pos, {x: this.pos.x + this.vel.x, y: this.pos.y + this.vel.y}], "white", 3);

    }
}

class Med {
    constructor() {
        this.place();
        this.radius = 4;
        this.collected = false;
        this.color = "blue";
    }
    place() {
        // TODO: locating a cell position
        // select random ring except innermost
        var rri = Math.round(1 + (Math.random()*(n_rings - 3)));
        console.log(rri);
        var length_ring = network.rings[rri].length;
        var rci = Math.floor(Math.random()*length_ring);
        this.pos = network.rings[rri][rci].get_pos();
    }
    check_collected() {
        // check for overlap
        if (distance(this.pos, player.pos) < (this.radius + player.radius)) {
            return true;
        } else {
            return false;
        }
    }
    effect() {
        // shrink player instantly back to start size
        player.radius = player_size;
    }
    update() {
        if (this.check_collected() && !this.collected) {
            this.effect();
            this.collected = true;
        }
    }
    render() {
        if (!this.collected) {
            draw_circ(this.radius, this.pos, this.color);
        }
    }
}

class Cell {
    constructor(ring_ind, ind_in_ring) {
        this.ring_ind = ring_ind
        this.ind_in_ring = ind_in_ring;
        // neighbours (also cells)
        this.neighbours = [];
        this.visited = false;
        this.color = "orange";
    }
    get_neighbours() {
        var outer_n = true;
        var inner_n = true;
        var this_ring = network.rings[this.ring_ind]
        var this_l = this_ring.length;
        var neighbours = [];
        // no outer neighbours if last ring
        if (this.ring_ind == n_rings - 1) { outer_n = false; }
        // no inner neighbours if first ring
        if (this.ring_ind == 0) { inner_n = false; }
        // determine outer neighbours
        if (outer_n) {
            var outer_l = network.rings[this.ring_ind + 1].length;
            // same ind if holes not doubled
            if (this_l == outer_l) {
                var outer_n_c = new Cell(this.ring_ind + 1, this.ind_in_ring);
                neighbours.push(outer_n_c);
            }
            // else: doubled - two outer neighbours
            else {
                var outer_n_c1 = new Cell(this.ring_ind + 1, this.ind_in_ring*2);
                var outer_n_c2 = new Cell(this.ring_ind + 1, this.ind_in_ring*2 + 1);
                neighbours.push(outer_n_c1);
                neighbours.push(outer_n_c2);
            }
        }
        // determine inner neighbours
        if (inner_n) {
            var inner_l = network.rings[this.ring_ind - 1].length;
            // - same ind if holes not doubled
            if (this_l == inner_l) {
                var inner_n_c = new Cell(this.ring_ind - 1, this.ind_in_ring);
                neighbours.push(inner_n_c);
            }
            // else: halved - different ind
            else {
                var inner_n_c = new Cell(this.ring_ind - 1, Math.floor(this.ind_in_ring/2));
                neighbours.push(inner_n_c);
            }
        }
        // determine sidewards neighbours
        if (this_l == 2) {
            // only left
            if (this.ind_in_ring + 1 == this_l) {
                var s_n = new Cell(this.ring_ind, 0);
            } else {
                var s_n = new Cell(this.ring_ind, this.ind_in_ring + 1);
            }
            neighbours.push(s_n);
        } else if (this_l > 2) {
            // left
            if (this.ind_in_ring + 1 == this_l) {
                var s_n1 = new Cell(this.ring_ind, 0);
            } else {
                var s_n1 = new Cell(this.ring_ind, this.ind_in_ring + 1);
            }
            // right
            if (this.ind_in_ring - 1 < 0) {
                var s_n2 = new Cell(this.ring_ind, this_l - 1);
            } else {
                var s_n2 = new Cell(this.ring_ind, this.ind_in_ring - 1);
            }
            neighbours.push(s_n1);
            neighbours.push(s_n2);
        }
        return neighbours;
    }
    get_rad() {
        var rads = equidistant_rads(network.rings[this.ring_ind].length);
        return rads[this.ind_in_ring];
    }
    get_pos() {
        // transform into coordinate
        // get distance from ring_ind
        var dist = this.ring_ind*ring_spacing + 0.5*ring_spacing;
        var rad = this.get_rad();
        return get_exact_coord(rad, dist, false);
    }
    render() {
        draw_circ(5, this.get_pos(), this.color);
    }
}

class MazeNetworkRandom {
    constructor() {
        // certain number of cells in each ring
        this.rings = [];
        // keeping track if number of holes should get doubled
        this.current_holes_n = 0;
        this.n_cells = 0;
        this.init_cells();
        this.current_cell = this.rings[0][0];
        this.previous_cell = this.rings[0][0];
        this.visited_cells = []; // stores the occupied cells
        this.stack = [];
        this.neighbour_cells = [];
        this.lines = [];
    }
    init_cells() {
        // BIG TODO: adjust so that n_cells == n_walls, always!
        // exception first ring: 4 cells even though no walls
        var n_holes = this.get_number_of_holes_ring(0);
        var cells = [];
        for (let index2 = 0; index2 < n_holes; index2++) {
            var cell = new Cell(0, index2);
            cells.push(cell);
            this.n_cells++;
        }
        this.rings.push(cells);
        for (let index = 1; index < n_rings; index++) {
            var cells = [];
            // get number of holes and set whether doubled
            var n_holes = this.get_number_of_holes_ring(index - 1);
            for (let index2 = 0; index2 < n_holes; index2++) {
                var cell = new Cell(index, index2);
                cells.push(cell);
                this.n_cells++;
            }
            this.rings.push(cells);
        }
    }
    set_closest_cell(pos) {
        // for debugging
        var current_dist = Infinity;
        var current_closest_cell = this.rings[0][0];
        for (let index = 0; index < this.rings.length; index++) {
            const ring = this.rings[index];
            for (let index2 = 0; index2 < ring.length; index2++) {
                const cell_pos = ring[index2].get_pos();
                if (distance(cell_pos, pos) < current_dist) {
                    current_dist = distance(cell_pos, pos);
                    current_closest_cell = ring[index2];
                }
            }
        }
        this.previous_cell = this.current_cell; // TODO: make copy?
        this.current_cell = current_closest_cell;
    }
    get_number_of_holes_ring(ring_ind) {
        // TODO some problem here - too many holes detected sometimes
        // how many holes (length r) fit onto circumference?
        var r = new Ring(ring_ind);
        var c = circumference_from_diameter(r.diameter);
        var n = Math.floor(c/ring_spacing);
        // make sure the walls fit as well
        while (n*ring_spacing + n*ring_thickness >= c) {
            n--;
        }
        // set to current holes
        if (this.current_holes_n == 0) {
            this.current_holes_n = n;
        } else if (n >= this.current_holes_n*2) {
            this.current_holes_n = this.current_holes_n*2;
        }
        return this.current_holes_n;
    }
    render() {
        // debug: render visited cells
        // check if last line already appended
        if (this.lines.length > 0) {
            var lastline = this.lines[this.lines.length-1];
            if (lastline[1].x != this.current_cell.get_pos().x || 
                lastline[1].y != this.current_cell.get_pos().y) {
                    var pos1 = this.previous_cell.get_pos();
                    var pos2 = this.current_cell.get_pos();
                    this.lines.push([pos1, pos2]);
            }
        } else {
            var pos1 = this.previous_cell.get_pos();
            var pos2 = this.current_cell.get_pos();
            this.lines.push([pos1, pos2]);
        }
        for (let index = 0; index < this.visited_cells.length; index++) {
            const visited_cell = this.visited_cells[index];
            visited_cell.color = "red";
            visited_cell.render();
            // draw line to current cell as well
            // if (index == this.visited_cells.length - 1) {
            //     var pos1 = visited_cell.get_pos();
            //     var pos2 = this.current_cell.get_pos();
            //     draw_line([pos1, pos2], "black");
            // }
        }
        for (let index = 0; index < this.lines.length; index++) {
            draw_line(this.lines[index], "black");
        }
        // debug: render current neighbour cells (set in maze.dig function)
        for (let index = 0; index < this.neighbour_cells.length; index++) {
            const neighbour = this.neighbour_cells[index];
            neighbour.color = "white";
            neighbour.render();
        }
        this.current_cell.color = "orange";
        this.current_cell.render();
    }
}

class MazeRandom {
    constructor(network) {
        // get info from network
        this.network = network;
        // generate the rings
        this.rings = [];
        this.go_back = false; // debug
        this.set_rings();
        this.init_walls();
        this.dig_all();
    }
    set_rings() {
        // already sets the holes accordingly
        for (let index = 0; index < n_rings; index++) {
            // add ring to array
            var r = new Ring(index);
            var c = circumference_from_diameter(r.diameter);
            // get number of holes
            var n = this.network.rings[index].length;
            // set holes
            // var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
            // hcoords = [];
            // r.set_holes_hcoords(hcoords);
            // add to rings array
            this.rings.push(r);
        }
    }
    init_walls() {
        // TODO: subdivision should occur at higher level
        for (let index = 1; index < this.rings.length; index++) {
            var n_walls = this.network.rings[index].length;
            var wall_coords = equidistant_wall_coords(n_walls);
            for (let index2 = 0; index2 < wall_coords.length; index2++) {
                var wall = new Wall(index - 1, index, wall_coords[index2]);
                this.rings[index].walls.push(wall);
            }
        }
    }
    remove_wall(ring_ind, wall_ind) {
        this.rings[ring_ind].walls[wall_ind].dug = true;
    }
    exclude_visited_neighbours(neighbours) {
        // TODO figure out why this returns one too many in case of going back
        var reduced_neighbours = [];
        for (let index = 0; index < neighbours.length; index++) {
            const nb = neighbours[index];
            var valid = true;
            for (let index2 = 0; index2 < network.visited_cells.length; index2++) {
                const vc = network.visited_cells[index2];
                if (nb.ring_ind == vc.ring_ind && nb.ind_in_ring == vc.ind_in_ring) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                reduced_neighbours.push(nb);
            }
        }
        return reduced_neighbours;
    }
    set_neighbours() {
        // first randomly select neighbouring cell
        var nbs = this.network.current_cell.get_neighbours();
        // make sure visited cells are excluded
        var nbs_MINUS_v = this.exclude_visited_neighbours(nbs);
        // debug: set the neighbour cells as property to network
        this.network.neighbour_cells = nbs_MINUS_v;
        return nbs_MINUS_v;
    }
    set_current_cell_to_last_valid_pos() {
        // goes through stack
        // determines number of valid neighbours
        // sets current cell as soon as there are >0 valid neighbours
        // also removes from stack
        // pop last element from stack
        network.stack.pop();
        const visited_cell = network.stack[network.stack.length - 1];
        network.previous_cell = network.current_cell; // TODO: make copy?
        network.current_cell = visited_cell;
        // determine neighbours
        var n_n = this.set_neighbours().length;
        if (n_n > 0) {
            // TODO figure out why 1 must be removed
            // stop search
            this.go_back = false;
            this.dig();
        }
    }
    set_visited_cell(c) {
        // first special case: no cells visited (start of program)
        if (network.stack.length == 0) {
            // set the current cell as visited
            c.visited = true;
            network.visited_cells.push(c);
            network.stack.push(c);
        } else {
            var last_visited = network.stack[network.stack.length - 1];
            if (c.ring_ind != last_visited.ring_ind || c.ind_in_ring != last_visited.ind_in_ring) {
                // set the current cell as visited no matter what
                c.visited = true;
                network.visited_cells.push(c);
                network.stack.push(c);
            }
        }
    }
    dig() {
        // removes walls or creates holes
        // TODO: put in correct place (debug)
        if (this.go_back) {
            this.set_current_cell_to_last_valid_pos();
            return;
        }
        // set cell as visited if not already set
        var c = network.current_cell;
        this.set_visited_cell(c);
        var n = this.set_neighbours();
        if (n.length == 0) {
            // set current cell to last location with options
            // but still store this cell as visited in list
            this.go_back = true;
            return;
        }
        // select randomly from cells
        var r_i = Math.floor(Math.random() * Math.floor(n.length));
        var rand_i = n[r_i];
        // introduce abbreviations
        var rri = rand_i.ring_ind;
        var cri = network.current_cell.ring_ind;
        var rir = rand_i.ind_in_ring;
        var cir = network.current_cell.ind_in_ring;
        var n_cells = network.rings[cri].length;
        var n_walls = this.rings[cri].walls.length;
        // sidewards: remove wall
        if (rri == cri && rri > 0) {
            if (rir == 0 && cir == network.rings[cri].length - 1) {
                // CASE 1
                // new cell right of old cell but wrapped around
                var dir = 1;
                var rwi = 0;
            } else if (rir == network.rings[cri].length - 1 && cir == 0) {
                // CASE 2
                // new cell left of old cell but wrapped around
                var dir = -1;
                var rwi = 0;
            } else if (rir < cir) {
                // CASE 3
                // new cell left of old cell and not wrapped around
                var dir = -1;
                var rwi = rir + 1;
            } else if (rir > cir) {
                // CASE 4
                // new cell right of old cell and not wrapped around
                var dir = 1;
                var rwi = rir;
            }
            // handle case: n_cells == 2*n_walls
            if (n_cells > n_walls) { // doubled
                if (dir == -1) {
                    if (rir == n_cells - 1) {
                        rwi = 0;
                    } else if (rir%2 == 1) {
                        rwi = Math.ceil(rir/2);
                    } else {
                        rwi = -1;
                    }
                } else {
                    if (rir%2 == 0) {
                        rwi = Math.floor(rir/2);
                    } else {
                        rwi = -1;
                    }
                }
            }
            // only set as dug if wall is to be removed
            if (rwi >= 0) {
                this.rings[rri].walls[rwi].dug = true;
            }
        }
        // through layers: add holes
        if (rri != cri) {
            // identify which ring
            var ri = Math.min(rri, cri);
            // setup holes coords to choose from
            var ring = network.rings[ri + 1];
            var n = ring.length;
            var c = this.rings[ri].get_c();
            var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
            // identify ind in ring
            if (rri < cri) {
                var ir = cir;
            } else {
                var ir = rir;
            }
            this.rings[ri].add_hcoord_to_hcoords(hcoords[ir]);
        }
        // set current cell to new position
        network.previous_cell = network.current_cell; // TODO: make copy?
        network.current_cell = rand_i;
    }
    dig_all() {
        var all_dug = false;
        while (!all_dug) {
            this.dig();
            if (network.visited_cells.length == network.n_cells) {
                all_dug = true;
            }
        }
        // set the last (escape) hole
        var ring = network.rings[network.rings.length-1];
        var n = ring.length;
        var c = this.rings[this.rings.length-1].get_c();
        var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
        var ir = Math.floor(Math.random()*hcoords.length);
        console.log(ir);
        this.rings[this.rings.length-1].add_hcoord_to_hcoords(hcoords[ir]);
    }
    update() {

        // if (network.visited_cells.length == network.n_cells) {
        //     console.log("ALL DUG");
        //     // set the last (escape) hole
        //     var ring = network.rings[network.rings.length-1];
        //     var n = ring.length;
        //     var c = this.rings[this.rings.length-1].get_c();
        //     var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
        //     var ir = Math.random()*hcoords.length;
        //     this.rings[this.rings.length-1].add_hcoord_to_hcoords(hcoords[ir]);
        // } else {
        //     this.dig();
        // }
    }
    render() {
        // maze consists of walls, rings
        // currently, walls belong to rings
        for (var i=0; i<this.rings.length; i++) {
            this.rings[i].render();
        }
    }
}

class Wall {
    constructor(ind1, ind2, rad, invalid=false, ring1, ring2) {
        this.ind1 = ind1; // inner ring to which wall is connected
        this.ind2 = ind2; // this ring to which wall is connected
        this.ring1 = ring1;
        this.ring2 = ring2;
        this.rad = rad; // angle in radians
        // derived vars
        this.r1 = ((ind1+1)*2*ring_spacing - ring_thickness)/2;
        this.r2 = ((ind2+1)*2*ring_spacing - ring_thickness)/2;
        this.coord1 = get_exact_coord(this.rad, this.r1, false);
        this.coord2 = get_exact_coord(this.rad, this.r2, false);
        this.invalid = invalid; // for debugging
        this.moving = false;
        this.dug = false;
    }
    set_validity() {
        // get them coords
        var coords1 = this.ring1.holes_coords;
        var coords2 = this.ring2.holes_coords;
        // if any overlap with hole of either inner or outer ring --> function sets validity
        this.invalid = false;
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
            this.rad -= 0.001;
            this.rad %= Math.PI*2; // make sure value stays within bounds
            this.update_coords();
        }
    }
    render() {
        var color = "black"
        if (!this.dug) {
            draw_line([this.coord1, this.coord2], color, ring_thickness);
        }
    }
}

class Ring {
    constructor(ind) {
        // lower inds --> smaller rings
        this.ind = ind;
        this.diameter = (ind+1)*2*ring_spacing - ring_thickness;
        this.holes_coords = [{start: 0, end: Math.PI*2}]; // init without holes
        this.edge_ball_coords = [];
        this.walls = [];
        this.rotated = false; // flag to store whether holes are rotated
        this.overlaps = [];
        this.overlaps_inner = [];
        this.debug_holes = [];
        this.debug_holes_inner = [];
        this.holes_overlap = false;
    }
    get_c() {
        // circumference from diameter
        return Math.PI*this.diameter;
    }
    set_holes(coords) {
        var coords_container = circumference_coords_given_coords(coords, this.diameter/2, ring_spacing, ring_thickness);
        this.holes_coords = coords_container[0];
        this.edge_ball_coords = coords_container[1];
    }
    add_coord_to_hcoords(coord) {
        // first remove placeholder:
        if (this.holes_coords.length == 1) {
            if (this.holes_coords[0].start == 0 && this.holes_coords[0].end == Math.PI*2) {
                this.holes_coords = [];
            }
        }
        var c = this.diameter*Math.PI;
        // flip coord (not sure why necessary)
        // coord = flip_rad(coord);
        var new_hcoord = start_end_from_coord(coord, c, ring_spacing, ring_thickness, true);
        this.holes_coords.push(new_hcoord);
        this.edge_ball_coords.push(new_hcoord.start);
        this.edge_ball_coords.push(new_hcoord.end);
    }
    add_hcoord_to_hcoords(hcoord) {
        // first remove placeholder:
        if (this.holes_coords.length == 1) {
            if (this.holes_coords[0].start == 0 && this.holes_coords[0].end == Math.PI*2) {
                this.holes_coords = [];
            }
        }
        this.holes_coords.push(hcoord);
        console.log(hcoord.start);
        // convert to coords
        var start = get_exact_coord(hcoord.start, this.diameter/2, false);
        var end = get_exact_coord(hcoord.end, this.diameter/2, false);
        this.edge_ball_coords.push(start);
        this.edge_ball_coords.push(end);
    }
    set_holes_hcoords(hcoords) {
        this.holes_coords = [];
        this.edge_ball_coords = [];
        for (let index = 0; index < hcoords.length; index++) {
            this.edge_ball_coords.push(hcoords[index].start);
            this.edge_ball_coords.push(hcoords[index].end);
        }
        this.holes_coords.push({start: hcoords[hcoords.length-1].end, end: hcoords[0].start});
        for (let index = 0; index < hcoords.length - 1; index++) {
            this.holes_coords.push({start: hcoords[index].end, end: hcoords[index + 1].start});
        }
    }
    /*
    set_holes_hcoords2() {
        if (this.holes_coords.length == 0) {
            this.holes_coords.push({start: 0, end: Math.PI*2});
        } else {
            this.edge_ball_coords = [];
            for (let index = 0; index < this.holes_coords.length; index++) {
                this.edge_ball_coords.push(this.holes_coords[index].start);
                this.edge_ball_coords.push(this.holes_coords[index].end);
            }
            this.holes_coords.push({start: this.holes_coords[hcoords.length-1].end, end: this.holes_coords[0].start});
            for (let index = 0; index < hcoords.length - 1; index++) {
                this.holes_coords.push({start: hcoords[index].end, end: hcoords[index + 1].start});
            }
        }
    }
    */
    add_wall(rad) {
        var lower_level_ring = rings[this.ind - 1];
        this.walls.push(new Wall(this.ind, this.ind-1, rad, true, this, lower_level_ring));
    }
    add_walls(lower_level_ring) {
        // /*
        // repeat finding random wall radians position until no overlap with holes
        var coords1 = this.holes_coords;
        var coords2 = lower_level_ring.holes_coords;
        var wall_rads = random_wall_rads(Math.ceil((this.ind + 1)*1.3));
        // var wall_rads = random_wall_rads(1);
        var validities = get_wall_validity(coords1, coords2, wall_rads);
        for (let index = 0; index < wall_rads.length; index++) {
            this.walls.push(new Wall(this.ind, this.ind-1, wall_rads[index], validities[index], this, lower_level_ring));
        }
        // */
    }
    update() {
        for (let index = 0; index < this.walls.length; index++) {
            this.walls[index].update();
        }
        for (let index = 0; index < this.walls.length; index++) {
            if (this.walls[index].invalid) {
                this.walls[index].moving = true;
                this.walls[index].update();
            } else {
                this.walls[index].moving = false;
            }
        }
    }
    render() {
        // convert holes_coords to line segment coords
        var line_seg_coords = coords_from_holes(this.holes_coords);
        // draw line segments
        for (let index = 0; index < line_seg_coords.length; index++) {
            draw_circ_segment(this.diameter/2, {x: canvas.width/2, y: canvas.height/2}, "black", ring_thickness, line_seg_coords[index]);
        }
        // draw edges
        for (let index = 0; index < this.edge_ball_coords.length; index++) {
            draw_circ(ring_thickness/2, this.edge_ball_coords[index], "black");
        }
        // draw walls
        for (let index = 0; index < this.walls.length; index++) {
            this.walls[index].render();
        }
    }
}

// instantiate objects
reset(min_rings);

// overall update function
function update() {
    // run changes (update objects)
    player.update();
    med.update();
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
    maze.render();
    // draw player
    player.render();
    med.render();
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

    // S, W --> debug increase/decrease player size
    if (e.keyCode == 87) { player.radius += 20*size_grow; }
    if (e.keyCode == 83 && player.radius > 1 + size_grow) { player.radius -= 20*size_grow; }

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
}

// start game loop
update();
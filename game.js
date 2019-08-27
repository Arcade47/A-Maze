var player_start_size;
var player_speed;
var ring_thickness;
var ring_spacing;
var n_rings;
var n_simulation_steps;
var growth_speed;
var player;
var maze;
var whole_len = Math.min(canvas.width, canvas.height);
var absolute;
var alpha;   
var beta;    
var gamma;
var debugball;

function new_level(set_n_rings) {
    n_rings = set_n_rings;
    ring_spacing = whole_len/(n_rings*2.5);
    player_start_size = 0.125*ring_spacing;
    ring_thickness = 0.125*ring_spacing;
    player_speed = 0.075*ring_spacing;
    n_simulation_steps = 30;
    growth_speed = 0.00025*ring_spacing;

    /*
    player_start_size = 5;
    player_speed = 4;
    ring_thickness = 5;
    ring_spacing = 40;
    n_rings = set_n_rings;
    n_simulation_steps = 30;
    growth_speed = 0.01;
    */

    player = new Player();
    maze = new Maze();
    maze.dig();
    debugball = new DebugBall();
}

document.addEventListener('keyup', keyup);
document.addEventListener('keydown', keydown);
window.addEventListener("deviceorientation", device_rotation);

class Ball {
    constructor(pos={x: center_coord.x, y: center_coord.y}, radius=ring_thickness/2) {
        this.pos = pos;
        this.radius = radius;
        this.color = "red";
    }
    render() {
        draw_circle_filled(this.pos, this.radius, this.color);
    }
}

class Player extends Ball {
    constructor() {
        super();
        this.color = "green";
        this.radius = player_start_size;
        this.moveleft = false;
        this.moveright = false;
        this.moveup = false;
        this.movedown = false;
        this.contact_points = [];
        this.shrinking = false;
        this.shrinking_percentage = 0;
    }
    set_contact_points_walls(ring_ind, dist, outside) {
        if (ring_ind > 0 && !outside) {
            for (let index = 0; index < maze.rings[ring_ind - 1].walls.length; index++) {
                const wall = maze.rings[ring_ind - 1].walls[index];
                if (!wall.dug) {
                    // get coord from distance and radians
                    var coord = rad_to_coord(wall.rad, dist);
                    this.contact_points.push(new Ball(coord));
                }
            }
        }
    }
    set_contact_points_rings(ring_ind) {
        var dist = (ring_ind + 1)*ring_spacing;
        // get radians from player pos
        var rad = coord_to_rad(this.pos);
        var coord = rad_to_coord(rad, dist);
        var in_hole = true;
        for (let index = 0; index < maze.rings[ring_ind].segments.length; index++) {
            const segment = maze.rings[ring_ind].segments[index];
            if (rad_between_rads(rad, [segment.start, segment.end])) {
                in_hole = false;
                break;
            }
        }
        if (!in_hole) {
            this.contact_points.push(new Ball(coord));
        }
    }
    set_contact_points_holes(ring_ind) {
        var dist = (ring_ind + 1)*ring_spacing;
        for (let index = 0; index < maze.rings[ring_ind].segments.length; index++) {
            const segment = maze.rings[ring_ind].segments[index];
            var start_coord = rad_to_coord(segment.start, dist);
            var end_coord = rad_to_coord(segment.end, dist);
            this.contact_points.push(new Ball(start_coord));
            this.contact_points.push(new Ball(end_coord));
        }
    }
    set_contact_points() {

        // reset
        this.contact_points = [];

        // identify relevant ring
        var center_dist = distance(this.pos, center_coord);
        var ring_ind = Math.floor(center_dist/ring_spacing);
        var outside = false;
        if (ring_ind > maze.n_rings-1) {
            ring_ind = maze.n_rings-1;
            outside = true;
        }
        var ring_dist = (ring_ind + 1)*ring_spacing;

        // walls
        this.set_contact_points_walls(ring_ind, center_dist, outside);

        // rings
        this.set_contact_points_rings(ring_ind);
        if (ring_ind > 0) {
            this.set_contact_points_rings(ring_ind - 1);
        }

        // holes
        this.set_contact_points_holes(ring_ind);
        if (ring_ind > 0) {
            this.set_contact_points_holes(ring_ind - 1);
        }

    }
    collision_resolution() {
        // velocity
        this.vel = {x: this.pos.x - this.old_pos.x, y: this.pos.y - this.old_pos.y};
        // simulation steps
        var step_x = this.vel.x / n_simulation_steps;
        var step_y = this.vel.y / n_simulation_steps;

        for (let index = 0; index < n_simulation_steps; index++) {

            // step position forward
            this.pos.x += step_x;
            this.pos.y += step_y;

            // did collision occur?
            this.set_contact_points();

            for (let index = 0; index < this.contact_points.length; index++) {
                const cp = this.contact_points[index];
                if (overlap(this, cp)) {
                    // correct for overlap
                    var dist_vect = {x: cp.pos.x - this.pos.x, y: cp.pos.y - this.pos.y};
                    var len = distance(this.pos, cp.pos);
                    var unit_dist_vect = {x: dist_vect.x/len, y: dist_vect.y/len};
                    var ball_overlap = (this.radius + cp.radius) - len;
                    var corr_vect = {x: unit_dist_vect.x*ball_overlap, y: unit_dist_vect.y*ball_overlap};
                    this.pos.x -= corr_vect.x;
                    this.pos.y -= corr_vect.y;
                }
            }
            
        }
    }
    shrink() {
        if (this.radius > player_start_size) {
            this.radius -= 4*growth_speed;
            var current_shrink_distance = this.radius - player_start_size;
            this.shrinking_percentage = current_shrink_distance/this.shrink_distance;
            if (this.radius <= player_start_size) {
                this.shrinking = false;
            }
        }
    }
    check_if_win() {
        var center_dist = distance(this.pos, center_coord);
        var ring_ind = Math.floor(center_dist/ring_spacing);
        
        if (ring_ind > maze.n_rings) {
            var new_n_rings = Math.min(21, maze.n_rings + 1);
            new_level(new_n_rings);
        }
    }
    check_if_lose() {
        // are remaining collisions?
        this.set_contact_points();

        for (let index = 0; index < this.contact_points.length; index++) {
            const cp = this.contact_points[index];
            if (overlap(this, cp)) {
                // get size of overlap
                var len = distance(this.pos, cp.pos);
                var ball_overlap = (this.radius + cp.radius) - len;
                if (ball_overlap > 0.2) {
                    var new_n_rings = Math.max(3, maze.n_rings - 1);
                    new_level(new_n_rings);
                }
            }
        }
    }
    update() {
        // store old pos
        this.old_pos = {x: this.pos.x, y: this.pos.y};

        if (!this.shrinking) {
            // increase player size
            this.radius += growth_speed;

            // react to keypresses
            if (this.moveleft) {
                this.pos.x -= player_speed;
            }
            if (this.moveright) {
                this.pos.x += player_speed;
            }
            if (this.moveup) {
                this.pos.y -= player_speed;
            }
            if (this.movedown) {
                this.pos.y += player_speed;
            }
            // collisions
            this.collision_resolution();

            // check for win/lose
            this.check_if_lose();
            this.check_if_win();
            
        } else {
            this.shrink();
        }
    }
}

class Item extends Ball {
    constructor(pos) {
        super();
        this.radius = 0.2*ring_spacing;
        this.pos = pos;
        this.opacity = 1;
        this.set_color();
        this.triggered = false;
        this.active = true;
    }
    set_color() {
        this.color = "rgba(0, 0, 255, " + String(this.opacity) + ")";
    }
    collision() {
        if (overlap(player, this)) {
            this.effect();
        }
    }
    effect() {
        // shrink player until start size
        if (player.shrinking == false) {
            player.shrink_start_abs = player.radius;
            player.shrink_distance = player.shrink_start_abs - player_start_size;
        }
        player.shrinking = true;
        this.triggered = true;
    }
    render() {
        if (this.active) {
            super.render();
        }
    }
    update() {
        if (!player.shrinking && this.triggered) {
            this.active = false;
        }
        if (player.shrinking && this.triggered) {
            this.opacity = player.shrinking_percentage;
            this.set_color(this.opacity);
        }
        if (this.active) {
            this.collision();
        }

    }
}

class DebugBall extends Ball {
    constructor() {
        super();
        this.pos = rad_to_coord(0, 100);
        this.color = "red";
    }
    update(a, b, g) {
        this.pos = rad_to_coord(degrees_to_rad(g), 100);
    }
}

class Wall {
    constructor(ind, rad) {
        // get 2 endpoints
        this.coords = [];
        this.rad = rad;
        var coord1 = rad_to_coord(rad, this.get_center_distance(ind));
        var coord2 = rad_to_coord(rad, this.get_center_distance(ind + 1));
        this.coords.push(coord1);
        this.coords.push(coord2);
        this.dug = false;
    }
    get_center_distance(ind) {
        return (ind + 1)*ring_spacing;
    }
    render() {
        if (!this.dug) {
            draw_line(this.coords, ring_thickness);
        }
    }
}

class Ring {
    constructor(ind) {
        // TODO add walls
        this.ind = ind;
        this.segments = [{start: 0, end: 2*Math.PI}];
        this.radius = (ind + 1)*ring_spacing;
        this.walls = [];
        this.whole = true;
    }
    add_walls(n_rings) {
        if (this.ind < n_rings - 1) {
            var walls_rads = get_walls_rads(this.radius, ring_thickness, ring_spacing);
            for (let index = 0; index < walls_rads.length; index++) {
                const wall_rad = walls_rads[index];
                this.walls.push(new Wall(this.ind, wall_rad));
            }
        }
    }
    dig_hole(ind_in_ring) {

        var n_holes_max = this.walls.length;
        // case: no outer walls --> look at next inner ring
        if (this.walls.length == 0) {
            n_holes_max = maze.rings[this.ind - 1].walls.length;
        }

        // identify position of hole
        var first_center_rad = (1.5*Math.PI) + (2*Math.PI)/(n_holes_max*2);
        var stepsize = (2*Math.PI)/n_holes_max;
        var center_rad = wrap_around(first_center_rad + (ind_in_ring*stepsize));

        var rad_distance = rad_dist(this.radius, ring_spacing);
        var h_start_val = center_rad - (0.5*rad_distance)
        var h_end_val = center_rad + (0.5*rad_distance)
        var hole_rads = {hstart: h_start_val, hend: h_end_val};
        
        // identify segment to cut open
        if (this.whole) {
            var start_value = hole_rads.hend;
            var end_value = hole_rads.hstart;
            this.segments = [{start: start_value, end: end_value}];
            this.whole = false;
            return;
        }
        var new_segments = [];
        for (let index = 0; index < this.segments.length; index++) {
            const segment = this.segments[index];
            if (rad_between_rads(center_rad, [segment.start, segment.end])) {
                var new_seg1 = {start: segment.start, end: hole_rads.hstart};
                var new_seg2 = {start: hole_rads.hend, end: segment.end};
                new_segments.push(new_seg1);
                new_segments.push(new_seg2);
            } else {
                new_segments.push(segment);
            }
        }
        this.segments = new_segments;

    }
    remove_wall(wall_ind) {
        this.walls[wall_ind].dug = true;
    }
    render() {
        // circle segments
        for (let index = 0; index < this.segments.length; index++) {
            const segment = this.segments[index];
            draw_circle_segment(center_coord, this.radius, segment, this.thickness);
        }
        // walls
        for (let index = 0; index < this.walls.length; index++) {
            const wall = this.walls[index];
            wall.render();
        }
    }
}

class Maze {
    constructor() {
        this.n_rings = n_rings;
        this.n_cells = 0;
        this.init_n_items();
        this.init_rings();
        this.ring_spacing = ring_spacing;
        this.visited_cells = [];
        this.paths = [];
        this.main_path = [];
        this.dead_ends = [];
        this.best_path = [];
        this.items = [];
    }
    init_n_items() {
        // 1 = m*5 + n
        // 1.5 = m*6 + n
        // 2 = m*7 + n
        // 2.5 = m*8 + n
        // ...
        // m = 0.5
        // 2.5 + n = 1 --> n = -1.5
        // MATH!
        this.n_items = Math.max(0, 0.5*this.n_rings - 1.5);
    }
    init_rings() {
        this.rings = [];
        for (let index = 0; index < this.n_rings; index++) {
            var ring = new Ring(index);
            ring.add_walls(this.n_rings);
            this.n_cells += ring.walls.length;
            this.rings.push(ring);
        }
        this.n_cells += this.rings[0].walls.length;
    }
    init_items(paths) {
        for (let item_ind = 0; item_ind < this.n_items; item_ind++) {

            // make sure there are enough paths
            if (paths.length-1 < item_ind) {
                break;
            }
            const de = paths[item_ind];
            const lc = de[de.length-1];

            // skip if in innermost ring
            if (lc.space_ind == 0) {
                continue;
            }

            this.items.push(new Item(lc.pos));
        }
    }
    dig() {
        var current_cell = this.get_random_start_cell();

        // while there are still unvisited cells
        while (this.visited_cells.length < this.n_cells) {

            var current_neighbours = current_cell.get_neighbours();
            var unvisited_neighbours = this.get_unvisited_neighbours(current_neighbours);
            
            var inserted_dead_end = false;  

            while (unvisited_neighbours.length == 0) {

                // dead end! go back --> only insert once
                if (!inserted_dead_end) {
                    this.dead_ends.push([]);
                    for (let index = 0; index < this.main_path.length; index++) {
                        const cell = this.main_path[index];
                        this.dead_ends[this.dead_ends.length-1].push(cell);
                        inserted_dead_end = true;
                    }
                    this.dead_ends[this.dead_ends.length-1].push(current_cell);
                }

                this.main_path.pop();
                current_cell = this.main_path[this.main_path.length-1];
                current_neighbours = current_cell.get_neighbours();
                unvisited_neighbours = this.get_unvisited_neighbours(current_neighbours);
            }

            // select random neighbour
            var next_cell = this.get_random_neighbour(unvisited_neighbours);

            // actually dig
            this.dig_place(current_cell, next_cell);

            // add to main path
            this.main_path.push(next_cell);

            // refresh current cell
            current_cell = next_cell;

            // if cell is valid
            this.visited_cells.push(current_cell);
        }

        // add to paths
        this.paths.push(this.main_path);
        for (let index = 0; index < this.dead_ends.length; index++) {
            const dead_end = this.dead_ends[index];
            this.paths.push(dead_end);
        }

        // dig exit
        this.dig_exit();
    }
    dig_place(cell1, cell2) {
        // identify separating maze component and erase
        if (cell1.space_ind != cell2.space_ind) {
            // 1. ring segment
            // 1.1 identify ring
            var ring_ind = Math.min(cell1.space_ind, cell2.space_ind);
            // 1.2 identify if same number of cells
            var n_cells1 = cell1.get_n_cells_in_space(cell1.space_ind);
            var n_cells2 = cell1.get_n_cells_in_space(cell2.space_ind);
            if (n_cells1 == n_cells2) {
                // 1.2.1 same number of cells in space
                this.rings[ring_ind].dig_hole(cell1.ind_in_space);
            } else {
                // 1.2.2 not same number of cells in space
                var ind_in_ring = Math.max(cell1.ind_in_space, cell2.ind_in_space);
                this.rings[ring_ind].dig_hole(ind_in_ring);
            }
        } else {
            // 2. wall
            // only dig in space where there are walls!
            if (cell1.space_ind > 0) {
                var n_cells = cell1.get_n_cells_in_space(cell1.space_ind);
                var ring_ind = cell1.space_ind - 1;
                var min_ind = Math.min(cell1.ind_in_space, cell2.ind_in_space);
                var max_ind = Math.max(cell1.ind_in_space, cell2.ind_in_space);
                if (min_ind == 0 && max_ind == n_cells-1) {
                    // 2.1 wrap around
                    this.rings[ring_ind].remove_wall(0);
                } else {
                    // 2.2 no wrap around
                    this.rings[ring_ind].remove_wall(max_ind);
                }
            }
        }
    }
    dig_exit() {

        var potential_paths = [];

        // 1. identify potential paths
        for (let index = 0; index < this.paths.length; index++) {
            const path = this.paths[index];
            for (let index2 = 0; index2 < path.length; index2++) {
                const cell = path[index2];
                if (cell.space_ind == this.rings.length-1) {
                    potential_paths.push(path);
                    break;
                }
            }
        }

        // 2. identify longest path
        var current_max_len = 0;
        var current_max_ind = -1;
        var potential_paths_copy = [];

        for (let index = 0; index < potential_paths.length; index++) {
            var path = potential_paths[index];
            var last_cell = path[path.length-1];
            potential_paths_copy.push(get_object_list_copy(path));
            while (last_cell.space_ind != this.rings.length-1) {
                path.pop();
                last_cell = path[path.length-1];
            }
            if (path.length > current_max_len) {
                current_max_len = path.length;
                current_max_ind = index;
            }
        }

        var best_path = potential_paths[current_max_ind];
        var furthest_cell = best_path[best_path.length-1];
        this.best_path = best_path;

        potential_paths_copy.splice(current_max_ind, 1);
        var potential_paths_lengths = this.get_path_deviation_lengths(potential_paths_copy, best_path);
        potential_paths_copy = sort_list_to_another_list(potential_paths_copy, potential_paths_lengths);

        // dig exit
        this.rings[this.rings.length-1].dig_hole(furthest_cell.ind_in_space);

        this.init_items(potential_paths_copy);

    }
    get_path_deviation_lengths(paths, main_path) {
        // return an array of integers: number of cells
        // that don't match main path
        var output = [];

        for (let index = 0; index < main_path.length; index++) {
            const mpc = main_path[index];
            for (let index2 = 0; index2 < paths.length; index2++) {
                const p = paths[index2];
                for (let index3 = 0; index3 < p.length; index3++) {
                    const pc = p[index3];
                    if (pc.space_ind != mpc.space_ind || pc.ind_in_space != mpc.ind_in_space) {
                        output.push(p.length - index3 + 1);
                        break;
                    }
                }
            }
        }

        return output;

    }
    get_random_neighbour(neighbours) {
        var l = neighbours.length;
        var rand = Math.random();
        // make absolutely shure there's no index error
        while (rand == 1) {
            rand = Math.random();
        }
        return neighbours[Math.floor(rand*l)];
    }
    get_random_start_cell() {
        var cell = new Cell(0, 0);
        var n_cells = cell.get_n_cells_in_space(0);
        var rand = Math.random();
        // make absolutely shure there's no index error
        while (rand == 1) {
            rand = Math.random();
        }
        var rand_ind_in_space = Math.floor(rand*n_cells);
        return new Cell(0, rand_ind_in_space);
    }
    get_unvisited_neighbours(neighbours) {
        var valid_neighbours = [];
        for (let index = 0; index < neighbours.length; index++) {
            const n = neighbours[index];
            var visited = false;
            for (let index2 = 0; index2 < this.visited_cells.length; index2++) {
                const v = this.visited_cells[index2];
                if (n.space_ind == v.space_ind && n.ind_in_space == v.ind_in_space) {
                    visited = true;
                }
            }
            if (!visited) {
                valid_neighbours.push(n);
            }
        }
        return valid_neighbours;
    }
    render() {
        // draw rings
        for (let index = 0; index < this.rings.length; index++) {
            const ring = this.rings[index];
            ring.render();
        }
        // draw walls
        // included in rings
        // draw items
        for (let index = 0; index < this.items.length; index++) {
            const item = this.items[index];
            item.render();
        }
    }
    update() {
        for (let index = 0; index < this.items.length; index++) {
            const item = this.items[index];
            item.update();
        }
    }
}

class Cell {
    constructor(space_ind, ind_in_space) {
        this.space_ind = space_ind;
        this.ind_in_space = ind_in_space;
        this.set_pos();
    }
    get_neighbours() {
        // neighbours inner, outer, sidewards
        var neighbours = [];

        // 1. vertical neighbours
        // 1.1 get indices of spaces
        var inner_neighbour_ind = -1;
        var outer_neighbour_ind = -1;
        if (this.space_ind > 0) {
            var inner_neighbour_ind = this.space_ind - 1;
        }
        if (this.space_ind < maze.rings.length-1) {
            var outer_neighbour_ind = this.space_ind + 1;
        }
        // 1.2 get indices in space
        var current_n_cells = this.get_n_cells_in_space(this.space_ind);
        // 1.2.1 inner
        if (inner_neighbour_ind > -1) {
            // compare number of cells
            var inner_n_cells = this.get_n_cells_in_space(inner_neighbour_ind);
            if (inner_n_cells < current_n_cells) {
                // 1.2.1.1 get index in space, halved
                var inner_ind_in_space = Math.floor(this.ind_in_space/2);
                neighbours.push(new Cell(inner_neighbour_ind, inner_ind_in_space));
            } else {
                // 1.2.1.2 get index in space, same inde
                neighbours.push(new Cell(inner_neighbour_ind, this.ind_in_space));
            }
        }
        // 1.2.2 outer
        if (outer_neighbour_ind > -1) {
            // compare number of cells
            var outer_n_cells = this.get_n_cells_in_space(outer_neighbour_ind);
            if (current_n_cells < outer_n_cells) {
                // 1.2.1.1 get index in space, doubled --> two cells
                var outer_ind_in_space1 = this.ind_in_space*2;
                var outer_ind_in_space2 = this.ind_in_space*2 + 1;
                neighbours.push(new Cell(outer_neighbour_ind, outer_ind_in_space1));
                neighbours.push(new Cell(outer_neighbour_ind, outer_ind_in_space2));
            } else {
                // 1.2.1.2 same number
                neighbours.push(new Cell(outer_neighbour_ind, this.ind_in_space));
            }
        }

        // 1. horizontal neighbours
        var n_cells_this_space = this.get_n_cells_in_space(this.space_ind);
        // 1.1 left neighbour
        if (this.ind_in_space == 0) {
            var left_ind_in_space = n_cells_this_space-1;
            neighbours.push(new Cell(this.space_ind, left_ind_in_space));
        } else {
            neighbours.push(new Cell(this.space_ind, this.ind_in_space-1));
        }
        // 1.2 right neighbour
        if (this.ind_in_space == n_cells_this_space-1) {
            var right_ind_in_space = 0;
            neighbours.push(new Cell(this.space_ind, right_ind_in_space));
        } else {
            neighbours.push(new Cell(this.space_ind, this.ind_in_space+1));
        }
        
        return neighbours;
    }
    get_n_cells_in_space(space_ind) {
        if (space_ind == 0) {
            return maze.rings[0].walls.length;
        }
        return maze.rings[space_ind - 1].walls.length;
    }
    set_pos() {
        // 1. get radians position
        // 1.1 get number of possible positions
        var n_rads_max = this.get_n_cells_in_space(this.space_ind);
        // 1.2 derive the rad
        var rad = rad_from_ind(n_rads_max, this.ind_in_space);
        // 2. get distance to center
        var dist_from_center = this.space_ind*ring_spacing + 0.5*ring_spacing;
        // 3. derive coord from rad and distance
        this.pos = rad_to_coord(rad, dist_from_center);
        // this.pos = coord_from_rad_distance(center_coord, dist_from_center, rad);
    }
}

new_level(3);

function update() {
    // updating objects
    player.update();
    maze.update();
    // draw functions
    render();
    // keep animation going
    window.requestAnimationFrame(update);
}

function render() {
    // refresh
    refresh_canvas("lightblue");
    // draw objects
    maze.render();
    player.render();
    debug_draw_text(alpha, 100);
    debug_draw_text(beta, 150);
    debug_draw_text(gamma, 200);
    debugball.render();
}

function keydown(e) {

    if (e.code == "ArrowUp") {
        player.moveup = true;
    }
    if (e.code == "ArrowDown") {
        player.movedown = true;
    }
    if (e.code == "ArrowLeft") {
        player.moveleft = true;
    }
    if (e.code == "ArrowRight") {
        player.moveright = true;
    }
}

function keyup(e) {
    
    if (e.code == "ArrowUp") {
        player.moveup = false;
    }
    if (e.code == "ArrowDown") {
        player.movedown = false;
    }
    if (e.code == "ArrowLeft") {
        player.moveleft = false;
    }
    if (e.code == "ArrowRight") {
        player.moveright = false;
    }
}

function device_rotation(e) {
    absolute = Math.round(e.absolute);
    alpha    = Math.round(e.alpha);
    beta     = Math.round(e.beta);
    gamma    = Math.round(e.gamma);
    debugball.update(alpha, beta, gamma);
}

// start the updating loop
update();
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
var no_virt_holes = false;
var complete = false;

// add event listeners
document.addEventListener("keydown", keydown);
document.addEventListener("keyup", keyup);

// prepare values for vars
function add_values() {
    // global vars that are initialized before
    n_rings = 5;            // rings of maze
    ring_spacing = 20;      // spacing of maze rings in px
    ring_thickness = 4;    // diameter of maze walls in px
    n_meds = 5;             // how many items need to be collected before allowed to exit
    size_grow = 0.005;          // increase of player diameter over a frame
    player_size = 3;        // initial player diameter in px
    player_speed = 2;       // movement speed
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

class MazeNetwork {

    // TODO layers are lists of nodes, i.e. network --> replace layers with network
    constructor(n_layers=n_rings, max_open_nodes=max_holes) {
        // input class variables
        this.n_layers = n_layers;
        this.max_open_nodes = max_open_nodes;

        // setup network (random node positions)
        this.n_nodes_network = this.generate_number_holes();
        this.n_dead_ends_network = this.generate_dead_end_numbers();
        this.n_virtual_nodes_network = this.generate_virtual_numbers();

        // create first network layout and loop randomly until valid
        // valid means: no overlaps in node order
        this.create_network_layout();
        while (this.overlap_network()) {
            this.create_network_layout();
        }

    }
    create_network_layout() {
        this.layers = this.setup_layers_with_nodes();
        this.setup_dead_ends();
        this.add_virt_nodes();
        this.connect_nodes();
    }
    generate_number_holes() {

        // following schmematic: 1,2,2,4,4,8,8,16,16,...,1
        var num = 2;
        var num_seq = [1];
        while (num_seq.length < this.n_layers - 1) {
            num_seq.push(num);
            num_seq.push(num);
            if (num < 8) { num *= 2; } // clip to upper bound
        }
        // handle case if even number --> add last num
        if (num_seq.length < this.n_layers) {
            num_seq.push(num);
        }
        num_seq[num_seq.length-1] = 1;
        return num_seq;
    }
    generate_dead_end_numbers() {
    
        // 0,1,0,1,0,2,0,4,0,8,0,...
        var num = 1;
        var num_seq = [];
        while (num_seq.length < this.n_layers - 1) {
            num_seq.push(0);
            num_seq.push(num);
            if (num_seq.length >= 4) {
                // clip to upper bound
                if (num < 2) {
                    num *= 2;
                } else {
                    num *= 0;
                }
            }
        }
        // handle case if even number --> add last num
        if (num_seq.length < this.n_layers) {
            num_seq.push(0);
        }
        // last layer: no dead ends
        num_seq[num_seq.length-1] = 0;
        
        return num_seq;
    }
    generate_virtual_numbers() {
        var num = 1;
        // assuming n >= 2
        var num_seq = [0,0];
        while (num_seq.length < this.n_layers - 1) {
            num_seq.push(num);
            num_seq.push(num);
            // clip to upper bound
            if (num < 4) {
                num *= 2;
            }
        }
        // handle case if even number --> add last num
        if (num_seq.length < this.n_layers) {
            num_seq.push(num);
        }
        return num_seq;
    }
    setup_layers_with_nodes() {
        var layers = [];
        for (let index = 0; index < this.n_layers; index++) {
            layers.push([]);
            if (index == 0 || index == this.n_layers - 1) {
                // only one node
                layers[layers.length-1].push(new Node(index, 0, this));
            } else {
                var n_nodes = this.n_nodes_network[index] - this.n_virtual_nodes_network[index];
                for (let index2 = 0; index2 < n_nodes; index2++) {
                    layers[layers.length-1].push(new Node(index, index2, this));
                }
            }
        }
        return layers;
    }
    setup_dead_ends() {
        // set the dead ends, starting from last layer moving upwards
        var last_layer_size = this.max_open_nodes; // keeping track of n nodes in last layer
        for (let index = this.layers.length-1; index >= 0; index--) {
            const layer = this.layers[index];
            var n_dead_ends = this.n_dead_ends_network[index];
            // adjust if last_layer_size smaller than current layer.length
            if (last_layer_size < layer.length && n_dead_ends < last_layer_size) {
                n_dead_ends = layer.length - last_layer_size;
            }
            // different case second-last layer: all but one are dead ends
            if (index == this.layers.length-2) {
                n_dead_ends = layer.length - 1;
            }
            var shuffle_array = [];
            for (let i = 0; i < n_dead_ends; i++) {
                shuffle_array.push(true);
            }
            for (let i = 0; i < layer.length - n_dead_ends; i++) {
                shuffle_array.push(false);
            }
            shuffle_array = shuffle(shuffle_array);
            for (let i = 0; i < shuffle_array.length; i++) {
                this.layers[index][i].dead_end = shuffle_array[i];
            }
            last_layer_size = layer.length;
        }
    }
    add_virt_nodes_layer(index, index2) {
        // add virtual node in in next layer if dead end or already virtual
        if (this.layers[index][index2].dead_end || this.layers[index][index2].virtual) {
            var virt_node = new Node(index + 1, index2, this);
            virt_node.virtual = true;
            this.layers[index + 1].splice(index2, 0, virt_node);
            // refresh all following IDs in this layer
            for (let index3 = index2 + 1; index3 < this.layers[index + 1].length; index3++) {
                this.layers[index + 1][index3].id++;
            }
            this.layers[index + 1][index2].set_parent(this.layers[index][index2]);
        }
    }
    add_virt_nodes() {
        for (let index = 0; index < this.layers.length - 1; index++) {
            for (let index2 = 0; index2 < this.layers[index].length; index2++) {
                this.add_virt_nodes_layer(index, index2);
            }
        }
    }
    connect_nodes() {

        for (let index = this.layers.length-1; index > 0; index--) {
            // going backwards, stopping at second layer
            // setting the parent of each node
            var upper_passage_count = 0;
            for (let i = 0; i < this.layers[index - 1].length; i++) {
                const node = this.layers[index - 1][i];
                if (!node.dead_end && !node.virtual) { upper_passage_count++; }
            }
            // identify virtual nodes in this layer
            var nonvirtual_count = this.layers[index].length;
            var nonvirtual_ids = [];
            for (let i = 0; i < this.layers[index].length; i++) {
                const node = this.layers[index][i];
                if (node.virtual) {
                    nonvirtual_count--;
                } else {
                    nonvirtual_ids.push(node.id);
                }
            }
        
            if (nonvirtual_count == upper_passage_count) { // one child per parent
                var upper_id = 0;
                for (let i = 0; i < nonvirtual_ids.length; i++) {
                    const nonvirtual_id = nonvirtual_ids[i];
                    for (let j = upper_id; j < this.layers[index - 1].length; j++) {
                        if (!this.layers[index - 1][j].dead_end && !this.layers[index - 1][j].virtual) {
                            this.layers[index][nonvirtual_id].set_parent(this.layers[index - 1][j]);
                            upper_id = j + 1;
                            break;
                        }
                    }
                }
        
            }
            else { // multiple children per parent
        
                var n_children_exceed = nonvirtual_count - upper_passage_count;
                // pick randomly from among parents which have multiple children - store number
                for (var i=0; i<n_children_exceed; i++) {
                    var rand_id = Math.floor(Math.random()*this.layers[index - 1].length);
                    // make sure no dead end or virtual node is targetted
                    var count = 0;
                    while ((this.layers[index - 1][rand_id].dead_end || this.layers[index - 1][rand_id].virtual) && count < 1000) {
                        var rand_id = Math.floor(Math.random()*this.layers[index - 1].length);
                        count++;
                    }
                    if (count >= 1000) {
                        console.log('problem here   ' + count);
                    }
                    this.layers[index - 1][rand_id].n_children++;
                }
                var lower_id = 0;
        
                for (let i = 0; i < this.layers[index - 1].length; i++) {
                    const upper_node = this.layers[index - 1][i];
                    var n_c = upper_node.n_children + 1; // + 1 because only exceeded children are noted --> one at least
                    if (!upper_node.dead_end && !upper_node.virtual) {
                        for (let j = 0; j < n_c; j++) {
                            while (this.layers[index][lower_id].virtual) {
                                lower_id++;
                            }
                            this.layers[index][lower_id].set_parent(this.layers[index - 1][i]);
                            lower_id++;
                        }
                    }
                }
        
            }
        }
    }
    overlap_network() {
        // returns false when there is no overlap
        for (let index = 0; index < this.layers.length - 1; index++) {
            for (let index2 = 0; index2 < this.layers[index].length - 1; index2++) {
                const node = this.layers[index][index2];
                // get rightmost child
                var right_id = node.get_id_of_rightmost_child();
                for (let index3 = index2 + 1; index3 < this.layers[index].length; index3++) {
                    const node_right = this.layers[index][index3];
                    // get leftmost child
                    var left_id = node_right.get_id_of_leftmost_child();
                    if (left_id < right_id) {
                        return true;
                    }
                }
            }
        }
        // no overlap problems
        return false;
    }
    render(current_node = [-1, -1]) {
        // makes sure the algo for the logical layout of maze is correct
        function get_coord(layer_ind, node_id) {
            var y_num = 10 + layer_ind*y_increase;
            var x_num = 10 + node_id*x_increase;
            return {x: x_num, y: y_num};
        }
        var y_increase = canvas.height/(this.layers.length*1.5) - 5;
        var x_increase = canvas.width/(this.max_open_nodes*2) - 5;
        var current_y = 10;
        var current_x = 10;
        var start_x = current_x;
        for (let index = 0; index < this.layers.length; index++) {
            const l = this.layers[index];
            current_x = start_x;
            for (let index2 = 0; index2 < l.length; index2++) {
                const n = l[index2];
                var n_color = "white";
                if (n.dead_end) { n_color = "red"; }
                if (n.virtual) { n_color = "grey"; }
                if (n.ind == current_node[0] && n.id == current_node[1]) { n_color = "green"};
                draw_circ_outline_debug(5, {x: current_x, y: current_y}, "black", n_color);
                // connect to parent
                for (let index3 = 0; index3 < n.parent.length; index3++) {
                    const p = n.parent[index3];
                    draw_line_debug([get_coord(p.ind, p.id), {x: current_x, y: current_y}], "black")
                }
                // set next position
                current_x += x_increase;
            }
            current_y += y_increase;
        }
    }
}

class Node {
    // virtual unit representing the maze structure
    constructor(ind, id, network) {
        this.ind = ind;
        this.id = id;
        this.parent = [];
        this.children = [];
        this.dead_end = false;
        this.virtual = false;
        this.n_children = 0;
        this.coord = 0;
        this.network = network;
    }
    set_parent(parent) {
        this.parent = [parent];
        // automatically sets children accordingly
        // but first search if child already included
        var child_included = false;
        for (let index = 0; index < this.network.layers[parent.ind][parent.id].children.length; index++) {
            var child = this.network.layers[parent.ind][parent.id].children[index];
            if (child.ind == this.ind && child.id == this.id) {
                child_included = true;
                break;
            }
        }
        if (!child_included) {
            this.network.layers[parent.ind][parent.id].children.push(this);
            // important flag to store number of children
            this.network.layers[parent.ind][parent.id].n_children = this.network.layers[parent.ind][parent.id].children.length;
        }
    }
    get_id_of_leftmost_child() {
        var min = Infinity;
        for (let index = 0; index < this.children.length; index++) {
            var child = this.children[index];
            if (child.id < min) {
                min = child.id;
            }
        }
        return min;
    }
    get_id_of_rightmost_child() {
        var max = -Infinity;
        for (let index = 0; index < this.children.length; index++) {
            var child = this.children[index];
            if (child.id > max) {
                max = child.id;
            }
        }
        return max;
    }
}

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
    set_coll_balls_along_ring(player_rad, ring_ind, dist) {
        var hole_check = overlap_hole_region(player_rad, rings[ring_ind].holes_coords, true);
            if (!hole_check[0]) { // exclude cases for inner ring if player is in innermost ring or holes
                // convert dist, rad to coord
                this.collision_balls.push(get_exact_coord(player_rad, dist, false));
            } else {
                // convert rads to coords
                this.collision_balls.push(get_exact_coord(flip_rad(hole_check[1].start), dist, false));
                this.collision_balls.push(get_exact_coord(flip_rad(hole_check[1].end), dist, false));
            }
    }
    set_coll_balls_along_walls(player_dist_center) {
        for (let index = 0; index < rings[this.current_ring_ind].walls.length; index++) {
            const wall = rings[this.current_ring_ind].walls[index];
            // convert to coords
            this.collision_balls.push(get_exact_coord(wall.rad, player_dist_center, false));
        }
    }
    set_closest_collision_balls() {
        // for static collision detection
        this.collision_balls = []; // refill array
        // two balls always except ring_ind == 0
        // get two distances from center (two rings)
        var dists = [];
        dists.push(this.current_ring_ind*ring_spacing - 0.5*ring_thickness);
        dists.push((this.current_ring_ind + 1)*ring_spacing - 0.5*ring_thickness);
        // also later needed: player distance from center
        var player_dist_center = distance(this.pos, center_coord);
        // get exact coords of collision_balls based on player angle
        var player_rad = coord_to_rad(this.pos);
        // start appending the virtual collision balls
        for (let index = 0; index < dists.length; index++) {
            if (this.current_ring_ind > 0 && this.current_ring_ind < n_rings) { // set bounds
                // 1. set collision balls along ring
                this.set_coll_balls_along_ring(player_rad, this.current_ring_ind + index - 1, dists[index]);
                // 2. set collision balls along walls
                this.set_coll_balls_along_walls(player_dist_center);
            } else if (this.current_ring_ind == 0 && index > 0) { // do not check artifical center ring, therefore index > 0
                // do not check walls
                this.set_coll_balls_along_ring(player_rad, 0, dists[index]);
            } else if (this.current_ring_ind >= n_rings) {
                // here also: do not check walls
                this.set_coll_balls_along_ring(player_rad, n_rings-1, (n_rings)*ring_spacing - 0.5*ring_thickness);
            }
        }
    }
    circular_move(CCW) {
        // get radians
        var pos_in_rad = coord_to_rad(this.pos);
        // get distance
        var dist = distance(this.pos, center_coord);
        // add radians
        var new_rad = pos_in_rad;
        if (CCW) { new_rad -= 0.05; } // TODO standardize on ring diameter
        if (!CCW) { new_rad += 0.05; }
        // convert back to coord
        this.pos = get_exact_coord(new_rad, dist, false);
    }
    resolve_collisions(old_pos) {
        // first reset position to old pos
        this.pos = {x: old_pos.x, y: old_pos.y};
        // split the velocity according to number of computations
        var vel_step = {x: this.vel.x/collision_steps, y: this.vel.y/collision_steps};
        // run the collision check collision_steps times --> only as long as no obstacle is reached
        var obstacle = false;
        for (let i = 0; i < collision_steps; i++) {
            if (!obstacle) {
                // step forward the position
                this.pos.x += vel_step.x;
                this.pos.y += vel_step.y;

                // circles thru the collision balls
                for (let j = 0; j < this.collision_balls.length; j++) {
                    const cball = this.collision_balls[j];
                    // check if overlap --> collision needs to be resolved
                    if (distance(this.pos, cball) < (this.radius + ring_thickness/2)) {
                        // assuming player movement caused collision --> vel > 0
                        var overlap = (this.radius + ring_thickness/2) - distance(this.pos, cball);
                        // scale velocity vector with overlap and substract
                        var unit_vec_vel = unit_vector(this.vel);
                        var scaled_vel = {x: unit_vec_vel.x*overlap, y: unit_vec_vel.y*overlap};
                        this.pos.x -= scaled_vel.x;
                        this.pos.y -= scaled_vel.y;
                        // stop the collision computation steps
                        obstacle = true;
                    }
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
        this.resolve_collisions(old_pos);

        // update current position in terms of ring indices
        this.set_current_ring_ind();

        // update virtual collision objects
        this.set_closest_collision_balls();

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

    }
}

class Maze {
    // TODO replaces old rings array
    constructor(network) {
        this.network = network;
        this.rings = [];
        this.randomly_init_rings();
        this.set_coords_of_innermost_ring();   
    }
    randomly_init_rings() {
        for (let index = 0; index < this.network.layers.length; index++) {
            var new_ring = new Ring(index, Math.random(), this.network.layers[index].length);
            this.rings.push(new_ring);
        }
    }
    set_coords_of_innermost_ring() {
        // set maximum rotation if adjacent rings have same number of holes
        // special case innermost ring: rotation does not matter
        var coords = circumference_coords_simpler(this.network.layers[0].length, start_rot);
        this.rings[0].set_holes(coords);
    }
}

class Wall {
    constructor(ind1, ind2, rad, invalid=false, ring1, ring2) {
        this.ind1 = ind1; // inner ring to which wall is connected
        this.ind2 = ind2; // outer ring to which wall is connected
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
        if (this.invalid) {
            color = "black";
        }
        draw_line([this.coord1, this.coord2], color, ring_thickness);
    }
}

class Ring {
    constructor(ind) {
        // lower inds --> smaller rings
        this.ind = ind;
        this.diameter = (ind+1)*2*ring_spacing - ring_thickness;
        this.holes_coords = [];
        this.edge_ball_coords = [];
        this.walls = [];
        this.rotated = false; // flag to store whether holes are rotated
        this.overlaps = [];
        this.overlaps_inner = [];
        this.debug_holes = [];
        this.debug_holes_inner = [];
        this.holes_overlap = false;
    }
    set_holes(coords) {
        var coords_container = circumference_coords_given_coords(coords, this.diameter/2, ring_spacing, ring_thickness);
        this.holes_coords = coords_container[0];
        this.edge_ball_coords = coords_container[1];
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
        // find new random wall rads as long as wall in hole
        // while (walls_in_holes(coords1, coords2, wall_rads)) {
        //     wall_rads = random_wall_rads(this.ind + 1);
        // }
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
    set_hole_problems() {
        // fills array with angles
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
        // draw walls
        for (let index = 0; index < this.walls.length; index++) {
            this.walls[index].render();
        }

        // debug: problematic wall overlaps
        for (let index = 0; index < this.overlaps.length; index++) {
            const rad = this.overlaps[index];
            // get outer coord
            var outer_coord = get_exact_coord(rad, -300, false);
            draw_line([center_coord, outer_coord], "red", 2);
        }
        // debug: problematic wall overlaps
        for (let index = 0; index < this.overlaps_inner.length; index++) {
            const rad = this.overlaps_inner[index];
            // get outer coord
            var outer_coord = get_exact_coord(rad, -300, false);
            draw_line([center_coord, outer_coord], "blue", 2);
        }

        // problematic holes
        for (let index = 0; index < this.debug_holes.length; index++) {
            draw_circ_segment(this.diameter/2, {x: canvas.width/2, y: canvas.height/2}, "green", ring_thickness/1, this.debug_holes[index]);
        }
        // problematic holes
        for (let index = 0; index < this.debug_holes_inner.length; index++) {
            draw_circ_segment(this.diameter/2, {x: canvas.width/2, y: canvas.height/2}, "red", ring_thickness/1, this.debug_holes_inner[index]);
        }
    }
}

// instantiate objects
add_values();
var player = new Player();

// construct maze according to rules
function create_network_layout() {
    // add nodes to the network
    var n_nodes_network = generate_number_holes(n_rings);
    var n_dead_ends_network = generate_dead_end_numbers(n_rings);
    var n_virtual_nodes_network = generate_virtual_numbers(n_rings);
    for (let index = 0; index < n_rings; index++) {
        layers.push([]);
        if (index == 0 || index == n_rings - 1) {
            // only one node
            layers[layers.length-1].push(new Node(index, 0));
        } else {
            var n_nodes = n_nodes_network[index] - n_virtual_nodes_network[index];
            for (let index2 = 0; index2 < n_nodes; index2++) {
                layers[layers.length-1].push(new Node(index, index2));
            }
        }
    }
    // set the dead ends, starting from last layer moving upwards
    var last_layer_size = max_holes; // keeping track of n nodes in last layer
    for (let index = layers.length-1; index >= 0; index--) {
        const layer = layers[index];
        // var n_dead_ends = Math.floor(dead_end_frac*layer.length);
        var n_dead_ends = n_dead_ends_network[index];
        // if (index == layers.length-2) {}
        // adjust if last_layer_size smaller than current layer.length
        if (last_layer_size < layer.length && n_dead_ends < last_layer_size) {
            n_dead_ends = layer.length - last_layer_size;
        }
        // different case second-last layer: all but one are dead ends
        if (index == layers.length-2) {
            n_dead_ends = layer.length - 1;
        }
        var shuffle_array = [];
        for (let i = 0; i < n_dead_ends; i++) {
            shuffle_array.push(true);
        }
        for (let i = 0; i < layer.length - n_dead_ends; i++) {
            shuffle_array.push(false);
        }
        shuffle_array = shuffle(shuffle_array);
        for (let i = 0; i < shuffle_array.length; i++) {
            layers[index][i].dead_end = shuffle_array[i];
        }
        last_layer_size = layer.length;
    }
    // set virtual nodes and connect them already
    function add_virt_nodes(index, index2) {
        // add virtual node in in next layer if dead end or already virtual
        if (layers[index][index2].dead_end || layers[index][index2].virtual) {
            var virt_node = new Node(index + 1, index2);
            virt_node.virtual = true;
            layers[index + 1].splice(index2, 0, virt_node);
            // refresh all following IDs in this layer
            for (let index3 = index2 + 1; index3 < layers[index + 1].length; index3++) {
                layers[index + 1][index3].id++;
            }
            layers[index + 1][index2].set_parent(layers[index][index2]);
        }
    }
    
    for (let index = 0; index < layers.length - 1; index++) {
        for (let index2 = 0; index2 < layers[index].length; index2++) {
            add_virt_nodes(index, index2);
        }
    }
    
    // connect the nodes
    for (let index = layers.length-1; index > 0; index--) {
        // going backwards, stopping at second layer
        // setting the parent of each node
        var upper_passage_count = 0;
        for (let i = 0; i < layers[index - 1].length; i++) {
            const node = layers[index - 1][i];
            if (!node.dead_end && !node.virtual) { upper_passage_count++; }
        }
        var nonvirtual_count = layers[index].length;
        var nonvirtual_ids = [];
        for (let i = 0; i < layers[index].length; i++) {
            const node = layers[index][i];
            if (node.virtual) {
                nonvirtual_count--;
            } else {
                nonvirtual_ids.push(node.id);
            }
        }
    
        if (nonvirtual_count == upper_passage_count) { // one parent one child
            // /*
    
            var upper_id = 0;
            for (let i = 0; i < nonvirtual_ids.length; i++) {
                const nonvirtual_id = nonvirtual_ids[i];
                for (let j = upper_id; j < layers[index - 1].length; j++) {
                    if (!layers[index - 1][j].dead_end && !layers[index - 1][j].virtual) {
                        layers[index][nonvirtual_id].set_parent(layers[index - 1][j]);
                        upper_id = j + 1;
                        break;
                    }
                }
            }
    
            // */
    
        } else { // multiple children per parent
    
            // /*
    
            var n_children_exceed = nonvirtual_count - upper_passage_count;
            // pick randomly from among parents which have multiple children - store number
            for (var i=0; i<n_children_exceed; i++) {
                var rand_id = Math.floor(Math.random()*layers[index - 1].length);
                // make sure no dead end or virtual node is targetted
                var count = 0;
                while ((layers[index - 1][rand_id].dead_end || layers[index - 1][rand_id].virtual) && count < 1000) {
                    var rand_id = Math.floor(Math.random()*layers[index - 1].length);
                    count++;
                }
                if (count >= 1000) {
                    console.log('problem here   ' + count);
                }
                layers[index - 1][rand_id].n_children++;
            }
            var lower_id = 0;
    
            // */
    
            // /*
    
            for (let i = 0; i < layers[index - 1].length; i++) {
                const upper_node = layers[index - 1][i];
                var n_c = upper_node.n_children + 1; // + 1 because only exceeded children are noted --> one at least
                if (!upper_node.dead_end && !upper_node.virtual) {
                    for (let j = 0; j < n_c; j++) {
                        while (layers[index][lower_id].virtual) {
                            lower_id++;
                        }
                        layers[index][lower_id].set_parent(layers[index - 1][i]);
                        lower_id++;
                    }
                }
            }
    
            // */
    
        }
    }
}
function overlap_network(layers) {
    // returns false when there is no overlap
    for (let index = 0; index < layers.length - 1; index++) {
        for (let index2 = 0; index2 < layers[index].length - 1; index2++) {
            const node = layers[index][index2];
            // get rightmost child
            right_id = node.get_id_of_rightmost_child();
            for (let index3 = index2 + 1; index3 < layers[index].length; index3++) {
                const node_right = layers[index][index3];
                // get leftmost child
                left_id = node_right.get_id_of_leftmost_child();
                if (left_id < right_id) {
                    return true;
                }
            }
        }
    }
    // no overlap problems
    return false;
}
function overlap_holes(this_coords, new_coords, inner_layer_ind) {

    var output = false; // flags if overlaps detected

    var inner_holes = holes_from_coords(this_coords);
    var this_holes = holes_from_coords(new_coords);

    for (let index = 0; index < inner_holes.length; index++) {
        var inner_coord = inner_holes[index];
        for (let index2 = 0; index2 < this_holes.length; index2++) {
        // for (let index2 = 1; index2 < 3; index2++) {
            var this_coord = this_holes[index2];
            var problem = overlap_two_holes(this_coord, inner_coord);

            if (problem[0]) {
                output = true;
            }
        }
    }

    return output;
}
function overlap_children(node) {
    // 1. get rad of parent_node.start
    var parent_layer_ind = node.ind;
    var parent_id = node.id;
    var parent_hcoords = holes_from_coords(rings[parent_layer_ind].holes_coords);
    var rad = parent_hcoords[parent_id].start;
    // 2. get children bounds
    var h_coords_children_layer = holes_from_coords(rings[parent_layer_ind + 1].holes_coords);
    // 2.1 case: only one child --> set larger bounds
    if (node.children.length == 1) {
        var child_id = node.children[0].id;
        var n_nodes_child_layer = network.layers[parent_layer_ind + 1].length;
        var next_child_id = (child_id + 1)%n_nodes_child_layer; // handle wrap-around case
        var min = h_coords_children_layer[child_id].start;
        var max = h_coords_children_layer[next_child_id].start;
    }
    // 2.2 case: more than one child
    else {
        // get correct indices
        var min_id = node.children[0].id;
        var max_id = node.children[node.children.length-1].id;
        // then convert holes coords to start/end of actual holes
        var min = h_coords_children_layer[min_id].start;
        var max = h_coords_children_layer[max_id].end;
    }
    var bounds = {start: min, end: max};
    // 3. test if overlap
    var test = rad_between_bounds(rad, bounds);
    // debug
    if (!test) {
        rings[parent_layer_ind].overlaps = [rad];
        rings[parent_layer_ind + 1].debug_holes_inner = [bounds];
    }
    return test;
    // }
}
function resolve_overlapping_holes() {

    const layer = network.layers[hole_problem_level];
    start_rot = (start_rot + 0.005)%1; // TODO set rotation angle proportional to circumference
    var coords = circumference_coords_simpler(layer.length, start_rot);
    rings[hole_problem_level].set_holes(coords);
    var this_hcoords = rings[hole_problem_level - 1].holes_coords;
    var inner_hcoords = rings[hole_problem_level].holes_coords;
    if (!overlap_holes(this_hcoords, inner_hcoords, hole_problem_level - 1)) {
        // var relevant_id = 0;
        for (let index = 0; index < network.layers[hole_problem_level - 1].length; index++) {
            var parent_node = network.layers[hole_problem_level - 1][index];
            if (parent_node.virtual) {
                relevant_id = parent_node.id;
            }
            if (hole_problem_level < n_rings - 1 && overlap_children(parent_node)) { // not last ring
                // success
                start_rot = Math.random();
                // already add the walls of this layer
                set_walls_one_layer(hole_problem_level);
                remove_virtuals_debug_layer(hole_problem_level);
                hole_problem_level++;
                rings[parent_node.ind].overlaps = [];
                rings[parent_node.ind + 1].debug_holes_inner = [];
                break;
            }
        }
    }
}
function rotate_last_ring() {
    // rotate last ring so long until 
    // find the one non-virtual node
    for (let index = 0; index < layers[layers.length-2].length; index++) {
        const node = layers[layers.length-2][index];
        if (!node.virtual) {
            while (!overlap_children(node)) {
                start_rot = (start_rot + 0.005)%1; // TODO set rotation angle proportional to circumference
                var coords = circumference_coords_simpler(layers[layers.length-2].length, start_rot);
                rings[layers[layers.length-2]].set_holes(coords);
            }
            break;
        }
    }
    // var parent_node = layers[layers.length-1][index];
    // if (overlap_children(parent_node)) {
    //     // TODO identify why not detected
    //     start_rot = Math.random();
    //     hole_problem_level++;
    //     rings[parent_node.ind].overlaps = [];
    //     rings[parent_node.ind + 1].debug_holes_inner = [];
    //     break;
    // }
}
function set_coords_of_innermost_ring(layers) {
    // set maximum rotation if adjacent rings have same number of holes
    // special case innermost ring: rotation does not matter
    var coords = circumference_coords_simpler(layers[0].length, start_rot);
    rings[0].set_holes(coords);
}
function set_walls_one_layer(index) {
    var tol = (rings[index].diameter/2)*0.00005;
    for (let index2 = 0; index2 < network.layers[index - 1].length; index2++) {
        const node = network.layers[index - 1][index2];
        // get position (simpler algorithm)
        // get correct indices
        var min_id = node.children[0].id;
        var max_id = node.children[node.children.length-1].id;
        // then convert holes coords to start/end of actual holes
        var hcoords_parent = holes_from_coords(rings[index - 1].holes_coords);
        var h_coords_children_layer = holes_from_coords(rings[index].holes_coords);
        var min = h_coords_children_layer[min_id].start;
        var max = h_coords_children_layer[max_id].end;
        // assuming parent start is always between bounds of children
        rings[index].add_wall((min - tol)%(2*Math.PI));
    }
}
function set_walls(layers) {
    for (let index = 1; index < layers.length; index++) { // do not start with innermost ring --> no walls!

        // add some tolerance to avoid unpleasing small overlaps
        // TODO error accumulation - small but bot quite correct

        // var circumference = 2*Math.PI*(rings[index].diameter)/2;
        // var tol = circumference*0.00005;
        var tol = (rings[index].diameter/2)*0.00005;

        for (let index2 = 0; index2 < layers[index - 1].length; index2++) {
            const node = layers[index - 1][index2];
            // get position (simpler algorithm)
            // get correct indices
            var min_id = node.children[0].id;
            var max_id = node.children[node.children.length-1].id;
            // then convert holes coords to start/end of actual holes
            var hcoords_parent = holes_from_coords(rings[index - 1].holes_coords);
            var h_coords_children_layer = holes_from_coords(rings[index].holes_coords);
            var min = h_coords_children_layer[min_id].start;
            var max = h_coords_children_layer[max_id].end;
            // assuming parent start is always between bounds of children
            rings[index].add_wall((min - tol)%(2*Math.PI));
        }
    }
}
function remove_virtuals() {
    for (let index = 0; index < rings.length; index++) {
        // BIGTODO
        
        var non_virtual_holes = [];
        var virtual_holes = [];

        // convert to holes coords
        var hcoords = holes_from_coords(rings[index].holes_coords);

        // omit virtual nodes
        for (let index2 = 0; index2 < hcoords.length; index2++) {
            const node = network.layers[index][index2];
            if (!node.virtual) {
                non_virtual_holes.push(hcoords[index2]);
            } else {
                virtual_holes.push(hcoords[index2]);
            }
        }

        // debug
        // rings[index].debug_holes = non_virtual_holes;
        // rings[index].debug_holes_inner = virtual_holes;

        rings[index].set_holes_hcoords(non_virtual_holes);
    }
}
function remove_virtuals_layer(index) {

    var non_virtual_holes = [];
    var virtual_holes = [];

    // convert to holes coords
    var hcoords = holes_from_coords(rings[index].holes_coords);

    // omit virtual nodes
    for (let index2 = 0; index2 < hcoords.length; index2++) {
        const node = layers[index][index2];
        if (!node.virtual) {
            non_virtual_holes.push(hcoords[index2]);
        } else {
            virtual_holes.push(hcoords[index2]);
        }
    }

    rings[index].set_holes_hcoords(non_virtual_holes);
}
function remove_virtuals_debug() {
    for (let index = 0; index < rings.length; index++) {
        
        var non_virtual_holes = [];
        var virtual_holes = [];

        // convert to holes coords
        var hcoords = holes_from_coords(rings[index].holes_coords);

        // omit virtual nodes
        for (let index2 = 0; index2 < hcoords.length; index2++) {
            const node = layers[index][index2];
            if (!node.virtual) {
                non_virtual_holes.push(hcoords[index2]);
            } else {
                virtual_holes.push(hcoords[index2]);
            }
        }

        // debug
        rings[index].debug_holes = non_virtual_holes;
        rings[index].debug_holes_inner = virtual_holes;

    }
}
function remove_virtuals_debug_layer(index) {

    var non_virtual_holes = [];
    var virtual_holes = [];

    // convert to holes coords
    var hcoords = holes_from_coords(rings[index].holes_coords);

    // omit virtual nodes
    for (let index2 = 0; index2 < hcoords.length; index2++) {
        const node = network.layers[index][index2];
        if (!node.virtual) {
            non_virtual_holes.push(hcoords[index2]);
        } else {
            virtual_holes.push(hcoords[index2]);
        }
    }

    // debug
    rings[index].debug_holes = non_virtual_holes;
    rings[index].debug_holes_inner = virtual_holes;

}
function resolve_overlaps_last_layer(parent_node) {
    const layer = network.layers[hole_problem_level]; // children layer
    start_rot = (start_rot + 0.005)%1; // TODO set rotation angle proportional to circumference
    var coords = circumference_coords_simpler(layer.length, start_rot);
    rings[hole_problem_level].set_holes(coords);
    var this_hcoords = rings[hole_problem_level - 1].holes_coords;
    var inner_hcoords = rings[hole_problem_level].holes_coords;
    if (!overlap_holes(this_hcoords, inner_hcoords, hole_problem_level - 1)) {

        // var parent_node = layers[hole_problem_level - 1][index];

        if (overlap_children(parent_node)) {
            // success
            // already add the walls of this layer
            set_walls_one_layer(hole_problem_level);
            remove_virtuals_debug_layer(hole_problem_level);
            hole_problem_level++;
            rings[parent_node.ind].overlaps = [];
            rings[parent_node.ind + 1].debug_holes_inner = [];
        } else {
            // console.log('no overlap children');
            rings[parent_node.ind + 1].debug_holes_inner = [this_hcoords[parent_node.id]];
        }

    }
}

// TODO erase old way of creating network
/*
var layers = [];
create_network_layout();
// make sure order of nodes in network is such that there are no overlaps
var network_invalid = overlap_network(layers);
while (network_invalid) {
    layers = [];
    create_network_layout();
    network_invalid = overlap_network(layers);
}
// debug: draw network
display_network_debug(layers, max_holes);
// */
var network = new MazeNetwork();

// add rings
// TODO delete; functions in constructor of maze class
// for (let index = 0; index < network.layers.length; index++) {
//     var new_ring = new Ring(index, Math.random(), network.layers[index].length);
//     rings.push(new_ring);
// }
// set_coords_of_innermost_ring(network.layers);

// make sure holes do not overlap
/*
while (hole_problem_level < n_rings) {
    resolve_overlapping_holes();
}
// */

// set_walls(layers);
// remove_virtuals();


// add the walls according to easier algorithm
// set_walls(layers);

// finally remove the virtuals
// remove_virtuals();

// overall update function
function update() {
    // run changes (update objects)
    player.update();
    // for (let index = 0; index < rings.length; index++) {
    //     rings[index].update();
    // }

    // /*
    if (hole_problem_level < n_rings - 1) {
        resolve_overlapping_holes();
    } else if (hole_problem_level == n_rings - 1) {
        // get relevant node
        var rel_ind = 0;
        for (let index = 0; index < network.layers[hole_problem_level - 1].length; index++) {
            const node = network.layers[hole_problem_level - 1][index];
            if (!node.virtual && !node.dead_end) {
                rel_ind = index;
                break;
            }
        }
        // TODO: make a unique function for last layer
        resolve_overlaps_last_layer(network.layers[hole_problem_level - 1][rel_ind]);

    } else if (!complete) {
        remove_virtuals();
        complete = true;
    }
    // */

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
    // draw player
    player.render();
    // draw network
    network.render();
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
        update();
        for (let index = 0; index < rings.length; index++) {
            for (let index2 = 0; index2 < rings[index].walls.length; index2++) {
                rings[index].walls[index2].moving = true;
            }
        }

    }

    // S, W --> debug increase player size
    if (e.keyCode == 87) {
        player.radius += 20*size_grow;
    }
    // S, W --> debug increase player size
    if (e.keyCode == 83 && player.radius > 1 + size_grow) { 
        player.radius -= 20*size_grow;
    }
    // Q, E --> make circular movements
    if (e.keyCode == 81) { player.CCWmove = true; }
    if (e.keyCode == 69) { player.CWmove = true; }

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
        for (let index = 0; index < rings.length; index++) {
            for (let index2 = 0; index2 < rings[index].walls.length; index2++) {
                rings[index].walls[index2].moving = false;
            }
        }
    }

    // Q, E --> make circular movements
    if (e.keyCode == 81) { player.CCWmove = false; }
    if (e.keyCode == 69) { player.CWmove = false; }
}

// start game loop
update();
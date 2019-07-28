var center_coord = {x: canvas.width/2, y: canvas.height/2};

function deg_to_rad(deg) {
    return deg*(Math.PI/180);
}

function flip_rad(rad) {
    return (rad + Math.PI) % (2*Math.PI);
}

function coords_diff(coord1, coord2) {
    var x1 = coord1.x;
    var x2 = coord2.x;
    var y1 = coord1.y;
    var y2 = coord2.y;
    var diff_x = x2 - x1;
    var diff_y = y2 - y1;
    return {x: diff_x, y: diff_y};
}

function distance(pos1,pos2) {
    var diff = coords_diff(pos1, pos2);
    // Pythagoras' theorem
    return Math.sqrt(diff.x*diff.x + diff.y*diff.y);
}

function coord_to_rad(coord) {
    // relative to center!
    var diff = coords_diff(coord, center_coord);
    return Math.atan2(diff.y, diff.x);
}

// trigonometric functions in JS are radians-based
function get_exact_coord(angle, dist, in_deg) {
    if (in_deg) {
        // TODO: flip angle
        var rad = deg_to_rad(angle); // convert
    } else {
        // flip angle (TODO: figure out why necessary)
        var rad = flip_rad(angle);
    }
    var x_val = dist*Math.cos(rad) + canvas.width/2;
    var y_val = dist*Math.sin(rad) + canvas.height/2;
    return {x: x_val, y: y_val};
}

function unit_vector(vec) {
    // get length
    var length = distance(vec, {x: 0, y: 0});
    // divide by length if not zero
    if (length > 0) {
        return {x: vec.x/length, y: vec.y/length};
    } else {
        return {x: 0, y: 0};
    }
}

function shuffle(a) {
    var j, x, i;
    for (i = a.length - 1; i > 0; i--) {
        j = Math.floor(Math.random() * (i + 1));
        x = a[i];
        a[i] = a[j];
        a[j] = x;
    }
    return a;
}

function flip_coords(coords) {
    var output = [];
    for (let index = 0; index < coords.length; index++) {
        const c = coords[index];
        output.push({start: flip_rad(c.start), end: flip_rad(c.end)});
    }
    return output;
}

function holes_from_coords(coords, flip = false) {

    var output = [];

    if (flip) {
        coords = flip_coords(coords);
    }

    // all holes but last one
    for (let j = 0; j < coords.length-1; j++) {
        var hole_start = coords[j].end %(Math.PI*2);
        var hole_end = coords[j+1].start %(Math.PI*2);
        output.push({start: hole_start, end: hole_end});
    }

    // last hole
    var hole_start = coords[coords.length-1].end %(Math.PI*2);
    var hole_end = coords[0].start %(Math.PI*2);
    output.push({start: hole_start, end: hole_end});
    
    return output;

}

function coords_from_holes(holes, flip=false) {
    // init with wrap around coord
    var output = [{start: holes[holes.length-1].end, end: holes[0].start}];

    if (flip) {
        coords = flip_coords(coords);
    }

    // reverse process from holes_from_coords
    for (let index = 1; index < holes.length; index++) {
        output.push({start: holes[index - 1].end, end: holes[index].start});
    }

    return output;
}

function generate_number_holes(n) {

    /*
    // following schmematic: 1,2,2,4,4,4,8,8,8,8,...,1
    var num_seq = [];
    var num = 1;
    var i = 0;
    while (num_seq.length < n) {
    // for (let index = 0; index < n; index++) {
        for (let index2 = 0; index2 < i + 1; index2++) {
            num_seq.push(num);
        }
        num *= 2;
        i += 1;
    }

    // */

    // following schmematic: 1,2,2,4,4,8,8,16,16,...,1
    var num = 2;
    var num_seq = [1];
    while (num_seq.length < n - 1) {
        num_seq.push(num);
        num_seq.push(num);
        if (num < 8) { num *= 2; } // clip to upper bound
    }
    // handle case if even number --> add last num
    if (num_seq.length < n) {
        num_seq.push(num);
    }
    num_seq[num_seq.length-1] = 1;
    return num_seq;
}

function generate_dead_end_numbers(n) {

    /*

    // 1,2,2,4,4,4,8,8,8,8 -->
    // 1,0,2,0,0,4,0,0,0,8 ...
    var num_seq = [];
    var num = 1;
    var i = 0;
    while (num_seq.length < n) {
    // for (let index = 0; index < n; index++) {
        for (let index2 = 0; index2 < i; index2++) {
            num_seq.push(0);
        }
        if (num_seq.length < n) {
            num_seq.push(num);
        }
        num *= 2;
        i += 1;
    }

    num_seq[0] = 0;
    num_seq[num_seq.length-1] = 0;

    // */

    // 0,1,0,1,0,2,0,4,0,8,0,...
    var num = 1;
    var num_seq = [];
    while (num_seq.length < n - 1) {
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
    if (num_seq.length < n) {
        num_seq.push(0);
    }
    // last layer: no dead ends
    num_seq[num_seq.length-1] = 0;
    
    return num_seq;
}

function generate_virtual_numbers(n) {
    var num = 1;
    // assuming n >= 2
    var num_seq = [0,0];
    while (num_seq.length < n - 1) {
        num_seq.push(num);
        num_seq.push(num);
        // clip to upper bound
        if (num < 4) {
            num *= 2;
        }
    }
    // handle case if even number --> add last num
    if (num_seq.length < n) {
        num_seq.push(num);
    }
    return num_seq;
}

function get_simpler_coords(coords) {
    // i.e. only starts
    var output = [];
    for (let index = 0; index < coords.length; index++) {
        const c = coords[index];
        output.push(c.start);
    }
    return output;
}
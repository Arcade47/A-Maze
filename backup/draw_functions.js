// inputs: radius, position on canvas and context

canv_w = 500;
canv_h = 500;

// TODO implement scrolling

function clear_canvas() {
	
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	
}

function set_canvas_bg(color) {
	
	ctx.fillStyle = color;
	ctx.fillRect(0,0,canv_w,canv_h);
	
}

function draw_grid() {
	
	for (h = 0; h < canvas.height; h += canvas.height/n_rows) {
		coords = [
			{x: 0, y: h},
			{x: canvas.width, y: h}
		];
		draw_poly(coords, "white");
	}
	
	for (v = 0; v < canvas.width; v += canvas.width/n_cols) {
		coords = [
			{x: v, y: 0},
			{x: v, y: canvas.height}
		];
		draw_poly(coords, "white");
	}
}

function refresh_canvas() {
	clear_canvas();
	draw_grid();
}

function draw_timed(n_frames) {
	if (visible_counter <= n_frames) {
		refresh_canvas();
		// draw stuff...
		visible_counter++;
	} else {
		refresh_canvas();
	}
}

function draw_line(coords, color, thickness=1) {
	ctx.beginPath();
	ctx.lineWidth=String(thickness);
	ctx.strokeStyle=color;
	ctx.moveTo(coords[0].x, coords[0].y);
	ctx.lineTo(coords[1].x, coords[1].y);
	ctx.stroke();
	ctx.closePath();
}

function draw_line_debug(coords, color, thickness=1) {
	ctx_l.beginPath();
	ctx_l.lineWidth=String(thickness);
	ctx_l.strokeStyle=color;
	ctx_l.moveTo(coords[0].x, coords[0].y);
	ctx_l.lineTo(coords[1].x, coords[1].y);
	ctx_l.stroke();
	ctx_l.closePath();
}

function draw_path(coords, color) {
	ctx.beginPath();
	ctx.lineWidth="1";
	ctx.strokeStyle=color;
	ctx.moveTo(coords[0].x, coords[0].y);
	for (let index = 1; index < coords.length; index++) {
		ctx.lineTo(coords[index].x, coords[index].y);
	}
	ctx.stroke();
	ctx.closePath();
}

function draw_circ(r, pos, color) {
	ctx.beginPath();
	ctx.fillStyle=color;
	ctx.arc(pos.x, pos.y, r, 0, 2*Math.PI);
	ctx.fill();
	ctx.closePath();
}

function draw_circ_outline(r, pos, outlinecolor, fillcolor) {
	ctx.beginPath();
	ctx.fillStyle=fillcolor;
	ctx.strokeStyle=outlinecolor;
	ctx.arc(pos.x, pos.y, r, 0, 2*Math.PI);
	ctx.fill();
	ctx.stroke();
	ctx.closePath();
}

function draw_circ_outline_debug(r, pos, outlinecolor, fillcolor) {
	ctx_l.beginPath();
	ctx_l.fillStyle=fillcolor;
	ctx_l.strokeStyle=outlinecolor;
	ctx_l.arc(pos.x, pos.y, r, 0, 2*Math.PI);
	ctx_l.fill();
	ctx_l.stroke();
	ctx_l.closePath();
}

function draw_circ_outline_transp(r, pos, outlinecolor, thickness) {
	ctx.beginPath();
	ctx.strokeStyle=outlinecolor;
	ctx.lineWidth = thickness;
	ctx.arc(pos.x, pos.y, r, 0, 2*Math.PI);
	ctx.stroke();
	ctx.closePath();
}

function draw_circ_segments(r, pos, outlinecolor, thickness, segments) {
	ctx.beginPath();
	ctx.strokeStyle=outlinecolor;
	ctx.lineWidth = thickness;
	for (let index = 0; index < segments.length; index++) {
		const start = segments[index].start;
		const end = segments[index].end;
		ctx.arc(pos.x, pos.y, r, start, end);
	}
	ctx.stroke();
	ctx.closePath();
}

function draw_circ_segment(r, pos, outlinecolor, thickness, segment) {
	ctx.beginPath();
	ctx.strokeStyle=outlinecolor;
	ctx.lineWidth = thickness;
	const start = segment.start;
	const end = segment.end;
	ctx.arc(pos.x, pos.y, r, start, end);
	ctx.stroke();
	ctx.closePath();
}

function draw_poly(coords, color) {
	ctx.beginPath();
	ctx.moveTo(coords[0].x, coords[0].y);
	ctx.strokeStyle=color;
	ctx.fillStyle=color;
	for (i=1; i<coords.length; i++) {
		ctx.lineTo(coords[i].x, coords[i].y);
	}
	ctx.lineTo(coords[coords.length-1].x, coords[coords.length-1].y);
	ctx.closePath();
	ctx.stroke();
	ctx.fill();
}

function draw_speech_bubble(text, pos) {
	// get measures from text
	ctx.font = "15px Arial";
	var textwidth = ctx.measureText(text).width;
	var textheight = ctx.measureText(text).height; // does not seem to work correctly
	var width = textwidth + 20;
	var height = 35;
	// create coords
	coords = [];
	coords.push([pos.x + 10, 					pos.y - 10]); // fixed
	coords.push([pos.x + 20, 					pos.y - 30]); // fixed
	coords.push([pos.x + 25 - width/2,		 	pos.y - 30]);
	coords.push([pos.x + 25 - width/2,			pos.y - 30 - height]);
	coords.push([pos.x + 25 + width/2,			pos.y - 30 - height]);
	coords.push([pos.x + 25 + width/2,			pos.y - 30]);
	coords.push([pos.x + 30, 					pos.y - 30]); // fixed
	coords.push([pos.x + 10, 					pos.y - 10]); // closing - first coord
	// adding colors
	ctx.strokeStyle="black";
	ctx.fillStyle="white";
	ctx.lineWidth = 3;
	ctx.beginPath();
	// assuming edge of bubble starts 10 px to right and above pos
	for (let index = 0; index < coords.length; index++) {
		const coord = coords[index];
		ctx.lineTo(...coord); //moveTo
	}
	ctx.closePath();
	ctx.stroke();
	ctx.fill();
	draw_canvas_text_flex(text, {x: pos.x + 25, y: pos.y - 25 - 0.5*height}, "black", 15, align="center");
}

function draw_vertex(coord) {
	
	ctx.beginPath();
	ctx.lineWidth="1";
	ctx.strokeStyle="blue";
	ctx.arc(coord.x,coord.y,5,0,2*Math.PI);
	ctx.stroke();
	ctx.closePath();
	
}

function draw_rect(coord, w, h, color) {
	ctx.beginPath();
	ctx.fillStyle=color;
	ctx.rect(coord.x, coord.y, w, h);
	ctx.fill();
}

function draw_rect_outline(coord, w, h, strokecolor, fillcolor) {
	ctx.beginPath();
	ctx.strokeStyle=strokecolor;
	ctx.fillStyle=fillcolor;
	ctx.rect(coord.x, coord.y, w, h);
	ctx.stroke();
	ctx.fill();
}

function draw_debug_text(string) {
	
	ctx.font = "40px Arial";
	ctx.strokeStyle="red";
	ctx.fillStyle="red";
	ctx.strokeText(string, 10, 50);
	ctx.fillText(string, 10, 50); 
	
}

function draw_debug_text_line(string, line) {
	
	ctx.font = "40px Arial";
	ctx.strokeStyle="red";
	ctx.fillStyle="red";
	ctx.fillText(string, 10, line*50); 
	
}

function draw_canvas_text_flex(string, pos, color, size, align="center") {
	
	ctx.font = String(size)+"px Arial";
	ctx.fillStyle=color;
	ctx.textAlign = align; 
	ctx.fillText(string, pos.x, pos.y);
	
}

function draw_canvas_text(string, pos) {
	
	ctx.font = "40px Arial";
	ctx.fillStyle="red";
	ctx.fillText(string, pos.x, pos.y);
	
}

function highlight_cell(cell, color) {
	var w = canvas.width/n_cols;
	var h = canvas.height/n_rows;
	draw_rect({x: cell.col*w, y: cell.row*h}, w, h, color);
}

function draw_problematic_cells(sec, pc, pcc) {
	
	// pc: problematic_cells
	// pcc: counter
	
	// draw problematic cells if any
	if (pc.length > 0 && pcc < sec*60) {
		for (var i=0; i<pc.length; i++) {
			highlight_cell({col: pc[i].col, row: pc[i].row}, "red");
			// move forward time
			pcc++;
			problematic_cells_counter = pcc;
		}
	}
	// reset if no problems or counter too long
	else {
		// set_values({probcells: [], probcount: 0});
		problematic_cells = [];
		problematic_cells_counter = 0;
	}
	
}

function draw_finished_polys(selectcolor, inactivecolor, list_of_polys, active_poly_ind) {
	
	for (var i=0; i<list_of_polys.length; i++) {
		// draw selected poly orange
		if (i == active_poly_ind) { var color = selectcolor; }
		else { var color = inactivecolor; }
		draw_poly(list_of_polys[i], color);
	}
	
}

function draw_polys_gameplay(color, list_of_polys) {
	
	for (var i=0; i<list_of_polys.length; i++) {
		draw_poly(list_of_polys[i], color);
	}
	
}

function draw_grid_occupancies(grid, color1, color2, color3, color4) {
	
	for (var i=0; i<grid.length; i++) {
		for (var j=0; j<grid[i].length; j++) {
			// if (grid[i][j] == -1) {highlight_cell({col: i, row: j}, color1);}
			if (grid[i][j] == 1) {highlight_cell({col: i, row: j}, color2);}
			else if (grid[i][j] == 2) {highlight_cell({col: i, row: j}, color3);}
			else if (grid[i][j] == 3) {highlight_cell({col: i, row: j}, color4);}
		}
	}
	
}

function draw_active_cells(grid, color) {
	for (var i=0; i<grid.length; i++) {
		for (var j=0; j<grid[i].length; j++) {
			if (grid[i][j] == -1) {highlight_cell({col: i, row: j}, color);}
		}
	}
}

function testdraw_active_lines(color, active_cells_draw) {
	
	for (var i=0; i<active_cells_draw.length; i++) {
		highlight_cell(active_cells_draw[i], color);
	}
	
}

function testdraw_polys_cells(color, list_of_polys_cells) {
	
	for (var i=0; i<list_of_polys_cells.length; i++) {
		for (var j=0; j<list_of_polys_cells[i].length; j++) {
			highlight_cell({col: list_of_polys_cells[i][j].col, row: list_of_polys_cells[i][j].row}, color);
		}
	}
	
}

function testdraw_cell_occupations(linecolor, selectedcolor, inactivecolor, polys_grid_state, active_poly_ind) {
	
	for (var i=0; i<polys_grid_state.length; i++) {
		for (var j=0; j<polys_grid_state[i].length; j++) {
			// active line
			if (polys_grid_state[i][j] == -1) {
				highlight_cell({col: i, row: j}, linecolor);
			}
			if (polys_grid_state[i][j] > 0 && polys_grid_state[i][j] != active_poly_ind+1) {
				highlight_cell({col: i, row: j}, inactivecolor);
			}
			if (polys_grid_state[i][j] > 0 && polys_grid_state[i][j] == active_poly_ind+1) {
				highlight_cell({col: i, row: j}, selectedcolor);
			}
		}
	}
	
}

function draw_active_vertices(active_poly) {
	
	for (var i=0; i<active_poly.length; i++) {
		draw_vertex(active_poly[i]);
	}
	
}

function draw_active_lines(active_poly) {
	
	var lines = lines_from_coords(active_poly, false);
	for (var i=0; i<lines.length; i++) {
		draw_poly(lines[i], "blue");
	}
	
}

function draw_cursor(pic, pos) {
	
	ctx.drawImage(
		pic,
		0,
		0,
		80,
		80,
		pos.x-10,
		pos.y-10,
		20,
		20
	);
	
}

function draw_highscores(names, scores) {
	set_canvas_bg("black");
	draw_canvas_text_flex("HIGHSCORES", 	{x: canvas.width/2, y: 40}, "white", 20, align="center");
	var ypos = 100;
	for (let index = 0; index < scores.length; index++) {
		draw_canvas_text_flex(String(index+1)+".", 	{x: canvas.width/2 - 180, y: ypos}, "white", 20, align="left");
		draw_canvas_text_flex(names[index], 		{x: canvas.width/2 - 160, y: ypos}, "white", 20, align="left");
		draw_canvas_text_flex(scores[index]+" €", 	{x: canvas.width - 180, y: ypos}, "white", 20, align="right");
		ypos += 30;
	}
	draw_canvas_text_flex("Press Enter to try again", 	{x: canvas.width/2, y: canvas.height - 30}, "white", 20, align="center");
}

function display_network(layers, max_holes, current_node = [-1, -1]) {
	// makes sure the algo for the logical layout of maze is correct
	function get_coord(layer_ind, node_id) {
		var y_num = 10 + layer_ind*y_increase;
		var x_num = 10 + node_id*x_increase;
		return {x: x_num, y: y_num};
	}
    var y_increase = canvas.height/(layers.length*1.5) - 5;
    var x_increase = canvas.width/(max_holes*2) - 5;
    var current_y = 10;
	var current_x = 10;
	var start_x = current_x;
    for (let index = 0; index < layers.length; index++) {
		const l = layers[index];
		current_x = start_x;
        for (let index2 = 0; index2 < l.length; index2++) {
			const n = l[index2];
			var n_color = "white";
			if (n.dead_end) { n_color = "red"; }
			if (n.virtual) { n_color = "grey"; }
			if (n.ind == current_node[0] && n.id == current_node[1]) { n_color = "green"};
			draw_circ_outline(5, {x: current_x, y: current_y}, "black", n_color);
			// connect to parent
			for (let index3 = 0; index3 < n.parent.length; index3++) {
				const p = n.parent[index3];
				draw_line([get_coord(p.ind, p.id), {x: current_x, y: current_y}], "black")
			}
			// set next position
			current_x += x_increase;
        }
        current_y += y_increase;
	}
	
	draw_debug_text(overlap_network(layers));

}

function display_network_debug(layers, max_holes, current_node = [-1, -1]) {
	// makes sure the algo for the logical layout of maze is correct
	function get_coord(layer_ind, node_id) {
		var y_num = 10 + layer_ind*y_increase;
		var x_num = 10 + node_id*x_increase;
		return {x: x_num, y: y_num};
	}
    var y_increase = canvas.height/(layers.length*1.5) - 5;
    var x_increase = canvas.width/(max_holes*2) - 5;
    var current_y = 10;
	var current_x = 10;
	var start_x = current_x;
    for (let index = 0; index < layers.length; index++) {
		const l = layers[index];
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
	
	// draw_debug_text(overlap_network(layers));

}

// control canvas prints - setup
var canvas = document.getElementById("GameCanvas");
var canvas_l = document.getElementById("LayoutCanvas");
canvas.width = canv_w;
canvas.height = canv_h;
canvas_l.width = canv_w;
canvas_l.height = canv_h;
ctx = canvas.getContext("2d");
ctx_l = canvas_l.getContext("2d");
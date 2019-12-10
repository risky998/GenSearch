/*
declare let input: string;
declare let output: any;
copy(text: string): text
print(text: string): undefined
*/

let amin = input;
let asteroidmap = amin.split("\n").map(l => l.split(""));
let resmap = amin.split("\n").map(l => l.split(""));

function gcd_two_numbers(x, y) {
	if (typeof x !== "number" || typeof y !== "number") return false;
	x = Math.abs(x);
	y = Math.abs(y);
	while (y) {
		var t = y;
		y = x % y;
		x = t;
	}
	return x;
}

function checkLineOfSight(x1, y1, x2, y2) {
	if (asteroidmap[y2][x2] !== "#") return false;
	if (asteroidmap[y1][x1] !== "#") return false;
	if (y2 === y1 && x2 === x1) return false;
	let dy = y2 - y1; // 4
	let dx = x2 - x1; // 2
	if (dx === dy && dy === 0) {
	} else {
		let gcd = gcd_two_numbers(dy, dx);
		dy /= gcd;
		dx /= gcd;
	}
	let x = x1 + dx;
	let y = y1 + dy;
	while (asteroidmap[y] && asteroidmap[y][x]) {
		if (asteroidmap[y][x] === "#") {
			if (y === y2 && x === x2) {
				return true;
			}
			return false;
		}
		y += dy;
		x += dx;
	}
}

let maxCount = 0;
let maxx = 0;
let maxy = 0;
for (let y = 0; y < asteroidmap.length; y++) {
	for (let x = 0; x < asteroidmap[0].length; x++) {
		let count = 0;
		for (let y2 = 0; y2 < asteroidmap.length; y2++) {
			for (let x2 = 0; x2 < asteroidmap[0].length; x2++) {
				count += +checkLineOfSight(x, y, x2, y2);
			}
		}
		resmap[y][x] = count;
		if (count > maxCount) {
			maxCount = count;
			maxx = x;
			maxy = y;
		}
	}
}

console.log(maxCount, maxx, maxy);
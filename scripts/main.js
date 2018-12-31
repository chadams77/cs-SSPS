window.SSPS = function () {

	this.vpw = window.innerWidth;
	this.vph = window.innerHeight;

	this.psim = new xSSPS(
		512, // Particle count
		512 // Render size
	);

	this.rcanvas = this.psim.canvas;

	this.canvas = document.createElement('canvas');
	this.canvas.width = this.vpw;
	this.canvas.height = this.vph;
	this.canvas.style.position = 'fixed';
	this.canvas.style.left = '0%';
	this.canvas.style.top = '0%';
	this.canvas.style.width = '100%';
	this.canvas.style.height = '100%';
	this.ctx = this.canvas.getContext('2d');
	this.ctx.imageSmoothingEnabled = true;
	this.ctx.imageSmoothingQuality = 'high';
	document.body.appendChild(this.canvas);

	this.keys = {};
};

SSPS.prototype.updateRender = function(dt) {

	this.psim.updateRender(this.keys, dt);

	this.ctx.clearRect(0, 0, this.vpw, this.vph);
	this.ctx.drawImage(this.rcanvas, 0, 0, this.rcanvas.width, this.rcanvas.height, this.vpw*0.5 - this.vph*0.5, 0, this.vph, this.vph);

	const lines = [];

	const mkOpt = (lst, i) => {
		i = i % lst.length;
		let ret = '';
		for (let j=0; j<lst.length; j++) {
			if (j === i) {
				ret += '>>' + lst[j] + '<< ';
			}
			else {
				ret += ' ' + lst[j] + ' ';
			}
		}
		return ret;
	};

	lines.push('SSPS by Chadams - ' + Math.round(1 / dt) + ' fps - ' + this.psim.PCOUNT + ' particles @ ' + this.psim.RSIZE + 'x' + this.psim.RSIZE);
	lines.push('[W] - Forward, [S] - Reverse, [ARROWS] - Look');
	lines.push('[E] - Grab, [SPACE] - Shoot Particle');
	lines.push('[R] - Cycle Render Mode: ' + mkOpt([1, 2, 3], this.psim.renderMode));
	lines.push('[1] - Cycle Particle Density: ' + mkOpt(this.psim.sDensity, this.psim.sDensityI));
	lines.push('[2] - Cycle Particle Viscosity: ' + mkOpt(this.psim.sVisc, this.psim.sViscI));
	lines.push('[3] - Cycle Particle Mass: ' + mkOpt(this.psim.sMass, this.psim.sMassI));
	lines.push('[4] - Cycle Particle Incompressiveness: ' + mkOpt(this.psim.sIncomp, this.psim.sIncompI));
	lines.push('[ESC] - Reset Particles & Camera')

	const fs = 14;
	const x = this.vpw*0.5 - this.vph*0.5 + 20;
	let y = 20 + (fs * 0.9);
	this.ctx.textAlign = 'left';
	this.ctx.font = fs + 'px Courier New';
	this.ctx.fillStyle = '#DDF';
	this.ctx.strokeStyle = 'rgba(0, 0, 0, 0.5)';

	for (let i=0; i<lines.length; i++) {
		this.ctx.fillText('' + lines[i], x, y);
		this.ctx.strokeText('' + lines[i], x, y);

		y += fs + 3;
	}


};

SSPS.prototype.start = function () {

	let lTime = Date.timeStamp();
	let lDt = 1/60;
	this.running = true;
	this.time = 0;

	document.body.addEventListener('keydown', (e) => {
		e = e || window.event;
		this.keys[e.keyCode] = true;
	});

	document.body.addEventListener('keyup', (e) => {
		e = e || window.event;
		this.keys[e.keyCode] = false;
	});

	const tick = () => {

		if (!this.running) {
			return;
		}

		this.updateViewport();

		const cTime = Date.timeStamp();
		const dt = (Math.max(Math.min(cTime - lTime, 1/10), 1/240) || (1/60)) * 0.1 + lDt * 0.9;
		lDt = dt;
		lTime = cTime;

		this.time += dt;

		this.updateRender(dt);

		requestAnimationFrame(tick);

	};

	requestAnimationFrame(tick);

};

SSPS.prototype.stop = function () {

	this.running = false;

	document.title = 'SSPS - Stopped';

};

SSPS.prototype.updateViewport = function() {

	if (this.vpw !== window.innerWidth || this.vph !== window.innerHeight) {
		this.vpw = window.innerWidth;
		this.vph = window.innerHeight;
		this.canvas.width = this.vpw;
		this.canvas.height = this.vph;
	}

}

Date.timeStamp = function() {
    return new Date().getTime() / 1000.0;
};
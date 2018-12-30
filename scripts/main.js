window.SSPS = function () {

	this.vpw = window.innerWidth;
	this.vph = window.innerHeight;

	this.psim = new xSSPS(512);

	document.body.appendChild(this.psim.canvas);
	this.psim.canvas.style.position = 'fixed';
	this.psim.canvas.style.zIndex = '1';
	this.psim.canvas.style.top = '0%';
	this.psim.canvas.style.left = 'calc(50vw - 50vh)';
	this.psim.canvas.style.width = this.psim.canvas.style.width = '100vh';
};

SSPS.prototype.updateRender = function(dt) {

	document.title = 'SSPS - ' + Math.floor(1 / dt) + ' fps';

	this.psim.updateRender(dt);

};

SSPS.prototype.start = function () {

	let lTime = Date.timeStamp();
	let lDt = 1/60;
	this.running = true;
	this.time = 0;

	document.body.addEventListener('click', (e) => {
		this.gindex += 1;
	});

	const tick = () => {

		if (!this.running) {
			return;
		}

		this.updateViewport();

		const cTime = Date.timeStamp();
		const dt = Math.max(Math.min(cTime - lTime, 1/10), 1/60) * 0.5 + lDt * 0.5;
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
	}

}
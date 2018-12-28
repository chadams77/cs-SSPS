window.SSPS = function () {

	this.vpw = window.innerWidth;
	this.vph = window.innerHeight;

	this.camera = new THREE.PerspectiveCamera(120, this.vpw / this.vph, 0.1, 3500);
	this.camera.position.z = 10;
	this.camera.lookAt(0, 0, 0);
	this.camera.updateProjectionMatrix();

	this.scene = new THREE.Scene();
	this.scene.background = new THREE.Color( 0x000000 );

	this.light = [];
	for (let i=0; i<5; i++) {
		this.light[i] = new THREE.PointLight( 0xffffff, 1, 1000 );
		this.light[i].position.set( 10000, 10000, 10000 );
		this.light[i].visible = false;
		this.scene.add( this.light[i] );
	}

	this.renderer = new THREE.WebGLRenderer({ alpha: true, antialias: false });
	this.renderer.setPixelRatio( window.devicePixelRatio );
	this.renderer.setSize( this.vpw, this.vph );
	document.body.appendChild(this.renderer.domElement);

	this.running = false;
	this.time = 0;

	this.psim = new xSSPS(this.scene, this.camera);
	this.gp = new THREE.Vector3(0, 0, 0);

	for (let i=0; i<2000; i++) {
		this.psim.add();
	}
};

SSPS.prototype.updateRender = function(dt) {

	document.title = 'SSPS - ' + Math.floor(1 / dt) + ' fps';

	const {grav, hl} = this.psim.updateRender(dt);

	this.gp.x += (grav[0].p.x - this.gp.x) * dt * 2;
	this.gp.y += (grav[0].p.y - this.gp.y) * dt * 2;
	this.gp.z += (grav[0].p.z - this.gp.z) * dt * 2;

	const ra = this.time/16 * Math.PI * 2;
	const ra2 = this.time/24 * Math.PI * 2;
	const rr = 5;
	this.camera.position.y = this.gp.y + Math.sin(ra2) * rr;
	this.camera.position.z = this.gp.z + Math.cos(ra) * Math.cos(ra2) * rr;
	this.camera.position.x = this.gp.x + Math.sin(ra) * Math.cos(ra2) * rr;
	this.camera.lookAt(this.gp.x, this.gp.y, this.gp.z);
	this.camera.updateProjectionMatrix();

	this.camera.up.y = this.gp.y;
	this.camera.up.z = this.gp.z + Math.cos(ra + Math.PI / 2) * rr;
	this.camera.up.x = this.gp.x + Math.sin(ra + Math.PI / 2) * rr;

	for (let i=0; hl && i<hl.length; i++) {
		this.light[i].position.set(hl[i].p.x, hl[i].p.y, hl[i].p.z);
		this.light[i].intensity = hl[i].heat * 0.005;
		this.light[i].visible = true;
	}
	for (let i=hl ? hl.length : 0; i<this.light.length; i++) {
		this.light[i].visible = false;
	}

};

SSPS.prototype.start = function () {

	let lTime = Date.timeStamp();
	let lDt = 1/60;
	this.running = true;
	this.time = 0;

	const tick = () => {

		if (!this.running) {
			return;
		}

		this.updateViewport();
		this.renderer.render( this.scene, this.camera );

		const cTime = Date.timeStamp();
		const dt = Math.max(Math.min(cTime - lTime, 1/24), 1/60) * 0.5 + lDt * 0.5;
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
		this.camera.aspect = this.vpw / this.vph;
		this.camera.updateProjectionMatrix();
		this.renderer.setSize( this.vpw, this.vph );
	}

}
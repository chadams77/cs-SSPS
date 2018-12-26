window.SSPS = function () {

	this.vpw = window.innerWidth;
	this.vph = window.innerHeight;

	this.camera = new THREE.PerspectiveCamera( 100, this.vpw / this.vph, 5, 3500 );

	this.scene = new THREE.Scene();
	this.scene.background = new THREE.Color( 0x000000 );

	this.renderer = new THREE.WebGLRenderer();
	this.renderer.setPixelRatio( window.devicePixelRatio );
	this.renderer.setSize( this.vpw, this.vph );

};

SSPS.prototype.updateRender = function(dt) {

	document.title = 'SSPS - ' + Math.floor(1 / dt) + ' fps';

};

SSPS.prototype.start = function () {

	let lTime = Date.timeStamp();

	const tick = () => {

		this.updateViewport();
		this.renderer.render( this.scene, this.camera );

		const cTime = Date.timeStamp();
		const dt = cTime - lTime;
		lTime = cTime;

		this.updateRender(dt);

		requestAnimationFrame(tick);

	};

	requestAnimationFrame(tick);

};

SSPS.prototype.updateViewport = function() {

	if (this.vpw !== window.innerWidth || this.vph !== window.innerHeight) {
		this.vpw = window.innerWidth;
		this.vph = window.innerHeight;
		camera.aspect = this.vpw / this.vph;
		camera.updateProjectionMatrix();
		renderer.setSize( this.vpw, this.vph );
	}

}
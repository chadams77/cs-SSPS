window.xSSPS = function(scene, camera) {

	const mkMaterial = (r,g,b,a) => {
		return new THREE.ShaderMaterial( {
			uniforms: {
				color: { value: new THREE.Vector4(r, g, b, a) },
			},
			vertexShader: document.getElementById( 'vertexShader1' ).textContent,
			fragmentShader: document.getElementById( 'fragmentShader1' ).textContent,
			transparent: true
		});
	};

	this.list = [];
	this.pgeom2 = new THREE.IcosahedronGeometry(1, 3);
	this.pgeom = new THREE.PlaneBufferGeometry(3, 3, 1, 1);
	//this.pgeom.computeFlatVertexNormals();
	this.mat = mkMaterial(1, 1, .1, 0.3); //new THREE.MeshStandardMaterial({ emissive: 0xFFFF10, color: 0xFFFF10, transparent: true });
	this.mat2 = mkMaterial(.1, .1, 1, 0.3); //new THREE.MeshStandardMaterial({ emissive: 0xFF0000, color: 0xFF0000, transparent: true });
	this.mat3 = new THREE.MeshStandardMaterial({ emissive: 0x000000, color: 0x000000, transparent: true });
	//this.mat.opacity = 0.1;
	this.mat.depthTest = true;
	this.mat.depthWrite = false;
	this.mat.blending = THREE.AdditiveBlending;
	this.mat.needsUpdate = true;
	this.mat2.fladShading = true;
	//this.mat2.opacity = 0.2;
	this.mat2.depthTest = true;
	this.mat2.depthWrite = false;
	this.mat2.needsUpdate = true;
	this.mat2.blending = THREE.AdditiveBlending;
	this.mat3.opacity = 1.0;
	this.mat3.depthTest = true;
	this.mat3.depthWrite = true;
	//this.mat3.blending = THREE.AdditiveBlending;
	this.mat3.needsUpdate = true;
	this.hSize = 0.25;
	this.mhSize = 2.5;
	this.nGravity = 15;
	this.gConst = 0.15;
	this.lFactor = 20.0;
	this.minMassGravity = 10;
	this.scene = scene;
	this.camera = camera;
	this.cmesh = [];
	for (let i=0; i<this.nGravity; i++) {
		const mesh = new THREE.Mesh(this.pgeom2, this.mat3);
		//this.scene.add(mesh);
		this.cmesh.push(mesh);
	}

};

const sphereDist = function(ax, ay, az, bx, by, bz, ar, br) {

	const dx = ax - bx, dy = ay - by, dz = az - bz;
	return Math.sqrt(dx*dx+dy*dy+dz*dz) - (ar+br);

};

xSSPS.prototype.updateRender = function(dt) {

	const hash3D = (x, y, z) => {
		const ix = (Math.floor(x/this.hSize)) + 1e5,
    	  	  iy = (Math.floor(y/this.hSize)) + 1e5,
    	  	  iz = (Math.floor(z/this.hSize)) + 1e5;
  		return ((ix * 10223) + (iy * 12919) + (iz * 16127)) % 1e6;
	};

	const mHash3D = (x, y, z) => {
		const ix = (Math.floor(x/this.mhSize)) + 1e5,
    	  	  iy = (Math.floor(y/this.mhSize)) + 1e5,
    	  	  iz = (Math.floor(z/this.mhSize)) + 1e5;
  		return ((ix * 10223) + (iy * 12919) + (iz * 16127)) % 1e6;
	};

	const hash = new Map();
	const mhash = new Map();

	// update position & angle, mark in hash
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		P.pos.x += P.vel.x * dt;
		P.pos.y += P.vel.y * dt;
		P.pos.z += P.vel.z * dt;
		P.ang.x += P.avel.x * dt;
		P.ang.y += P.avel.y * dt;
		P.ang.z += P.avel.z * dt;
		P.hkey = hash3D(P.pos.x, P.pos.y, P.pos.z);
		P.mhkey = mHash3D(P.pos.x, P.pos.y, P.pos.z);

		// insert into collision hash
		if (!hash.has(P.hkey)) {
			hash.set(P.hkey, [P]);
		}
		else {
			hash.get(P.hkey).push(P);
		}

		// insert into gravity hash
		if (!mhash.has(P.mhkey)) {
			mhash.set(P.mhkey, [P]);
		}
		else {
			mhash.get(P.mhkey).push(P);
		}

		// update THREE.js mesh
		P.mesh.position.set(P.pos.x, P.pos.y, P.pos.z);
		//P.mesh.rotation.set(P.ang.x, P.ang.y, P.ang.z, 'XYZ');
		P.mesh.quaternion.copy(this.camera.quaternion);
		P.mesh.scale.set(P.lr, P.lr, P.lr);
		P.mesh.updateMatrix(true);
	}

	// compute gravity hash
	const gravity = [];
	for (let [hk, list] of mhash) {
		const gp = new THREE.Vector3(0, 0, 0);
		let tmass = 0.0;
		for (let i=0; i<list.length; i++) {
			const P = list[i];
			const mass = P.mass;
			gp.x += P.pos.x * mass;
			gp.y += P.pos.y * mass;
			gp.z += P.pos.z * mass;
			tmass += mass;
		}
		if (tmass > this.minMassGravity) {
			gp.x /= tmass;
			gp.y /= tmass;
			gp.z /= tmass;
			gravity.push({
				p: gp,
				mass: tmass
			});
		}
	}
	gravity.sort((a, b) => (b.mass - a.mass));
	gravity.length = Math.min(gravity.length, this.nGravity);


	let gtmass = gravity[0].mass;
	let gtp = new THREE.Vector3(gravity[0].p.x, gravity[0].p.y, gravity[0].p.z);
	for (let i=1; i<gravity.length; i++) {
		gtp.x = (gtp.x * gtmass + gravity[i].p.x * gravity[i].mass) / (gravity[i].mass + gtmass);
		gtp.y = (gtp.y * gtmass + gravity[i].p.y * gravity[i].mass) / (gravity[i].mass + gtmass);
		gtp.z = (gtp.z * gtmass + gravity[i].p.z * gravity[i].mass) / (gravity[i].mass + gtmass);
		gtmass += gravity[i].mass;
	}

	if (!this.gtp) {
		this.gtp = gtp;
		this.gtmass = gtmass;
	}
	else {
		this.gtp.x += (gtp.x - this.gtp.x) * dt * 10;
		this.gtp.y += (gtp.y - this.gtp.y) * dt * 10;
		this.gtp.z += (gtp.z - this.gtp.z) * dt * 10;
		this.gtmass += (gtmass - this.gtmass) * dt * 6;
	}

	for (let i=0; i<this.cmesh.length; i++) {
		this.cmesh[i].scale.set(0.0001, 0.0001, 0.0001);
	}
	for (let i=0; i<1; i++) {
		this.cmesh[i].position.set(
			this.gtp.x, this.gtp.y, this.gtp.z
		);
		const gscale = 0.1 + Math.pow(this.gtmass, 1/3) / 10;
		this.cmesh[i].scale.set(
			gscale, gscale, gscale
		);
		this.cmesh[i].updateMatrix(true);
		this.cmesh[i].updateWorldMatrix(true);
	}

	// update velocities with gravity
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		for (let j=0; j<gravity.length; j++) {
			const G = gravity[j];
			const dx = G.p.x - P.pos.x,
				  dy = G.p.y - P.pos.y,
				  dz = G.p.z - P.pos.z;
			const dlenSq = (dx*dx + dy*dy + dz*dz);
			if (dlenSq > 0.01) {
				const dlen1 = Math.sqrt(dlenSq);
				const dlen2 = G.mass / (Math.max(dlenSq, 20) * 0.1);
				P.vel.x += dlen2 * (dx / dlen1) * dt * this.gConst;
				P.vel.y += dlen2 * (dy / dlen1) * dt * this.gConst;
				P.vel.z += dlen2 * (dz / dlen1) * dt * this.gConst;
			}
		}
	}

	// handle collisions
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		let ix, iy, iz;
		let dx, dy, dz, dlen;
		const tested = new Map();
		for (ix=-1; ix<=1; ix++) {
			for (iy=-1; iy<=1; iy++) {
				for (iz=-1; iz<=1; iz++) {
					const ihk = hash3D(
						P.pos.x + ix * P.lr,
						P.pos.y + iy * P.lr,
						P.pos.z + iz * P.lr
					);
					if (hash.has(ihk)) {
						const list2 = hash.get(ihk);
						for (let j=0; j<list2.length; j++) {
							const jP = list2[j];
							if (jP.id > P.id && !tested.has(jP.id)) {
								tested.set(jP.id, true);
								const sdist = sphereDist(
									P.pos.x, P.pos.y, P.pos.z,
									jP.pos.x, jP.pos.y, jP.pos.z,
									P.lr, jP.lr
								);
								if (sdist < 0) {
									dx = P.pos.x - jP.pos.x;
									dy = P.pos.y - jP.pos.y;
									dz = P.pos.z - jP.pos.z;
									dlen = Math.sqrt(dx*dx+dy*dy+dz*dz);
									if (dlen > 0.0001) {
										dx /= dlen;
										dy /= dlen;
										dz /= dlen;
										const damp = 0.5 * (P.ldamp + jP.ldamp);
										const force = this.lFactor * Math.pow(-sdist / Math.min(P.lr, jP.lr), 2.0) * (P.mass + jP.mass);

										P.vel.x += dx * force * dt / P.mass;
										P.vel.y += dy * force * dt / P.mass;
										P.vel.z += dz * force * dt / P.mass;
										jP.vel.x -= dx * force * dt / jP.mass;
										jP.vel.y -= dy * force * dt / jP.mass;
										jP.vel.z -= dz * force * dt / jP.mass;
										P.pos.x += dx * force * dt * dt / P.mass;
										P.pos.y += dy * force * dt * dt / P.mass;
										P.pos.z += dz * force * dt * dt / P.mass;
										jP.pos.x -= dx * force * dt * dt / jP.mass;
										jP.pos.y -= dy * force * dt * dt / jP.mass;
										jP.pos.z -= dz * force * dt * dt / jP.mass;

										const ivx = (P.vel.x * P.mass + jP.vel.x * jP.mass) / (P.mass + jP.mass);
										const ivy = (P.vel.y * P.mass + jP.vel.y * jP.mass) / (P.mass + jP.mass);
										const ivz = (P.vel.z * P.mass + jP.vel.z * jP.mass) / (P.mass + jP.mass);

										const dampaf = Math.pow(P.ldamp, dt);
										const dampbf = Math.pow(jP.ldamp, dt);

										P.vel.x = P.vel.x * dampaf + ivx * (1 - dampaf);
										P.vel.y = P.vel.y * dampaf + ivy * (1 - dampaf);
										P.vel.z = P.vel.z * dampaf + ivz * (1 - dampaf);
										jP.vel.x = jP.vel.x * dampaf + ivx * (1 - dampbf);
										jP.vel.y = jP.vel.y * dampaf + ivy * (1 - dampbf);
										jP.vel.z = jP.vel.z * dampaf + ivz * (1 - dampbf);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return [{p: this.gtp, mass: this.gtmass}];

};

xSSPS.prototype.add = function(inArgs) {

	const args = inArgs || {};

	const rat = Math.PI * Math._random(), rap = Math.PI * 2.0 * Math._random();
	const rr = Math._random() * 60;
	const id = 1e9 + Math.floor(Math._random() * 1e9);
	const pos = args.pos ? args.pos : new THREE.Vector3(Math.cos(rat) * Math.sin(rap) * rr, Math.sin(rat) * Math.sin(rap) * rr, Math.cos(rap) * rr);

	const type = (rr < 35) ? 1 : 0;

	const ang = args.ang || new THREE.Vector3(Math._random()*50, Math._random()*50, Math._random()*50);
	const mesh = new THREE.Mesh(this.pgeom, [this.mat, this.mat2, this.mat3][type]);

	mesh.position.set(pos.x, pos.y, pos.z);
	mesh.rotation.set(ang.x, ang.y, ang.z, 'XYZ');

	const obj = {
		id: id,
		group: null,
		pos: pos,
		ang: ang,
		vel: args.vel || new THREE.Vector3(0, 0, 0),
		avel: args.avel || new THREE.Vector3(0, 0, 0),
		mass: [0.25, 1.0, 8.0][type],
		mesh: mesh,
		hr: 0.05,
		lr: ([0.5, 0.5, 4.0][type]) * 0.5,
		hdamp: 0.1,
		ldamp: [0.75, 0.55, 0.05][type],
		hkey: 0
	};

	this.scene.add(mesh);

	this.list.push(obj);

};
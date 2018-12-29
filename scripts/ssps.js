window.xSSPS = function(scene, camera) {

	this.list = [];
	this.pgeom2 = new THREE.IcosahedronGeometry(0.5, 0);
	this.pgeom = new THREE.PlaneBufferGeometry(2, 2, 1, 1);
	this.hSize = 1.0;
	this.fLen = this.hSize;
	this.mhSize = 1024;
	this.mhMinSize = 1;
	this.restMass = 20;
	this.nGravity = 15;
	this.gConst = 0.001;
	this.lFactor = 20.0;
	this.minMassGravity = 10;
	this.minHeatLight = 10;
	this.nHeat = 5;
	this.scene = scene;
	this.camera = camera;

	this.seedTypes();

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

	const hash3Di = (x, y, z) => {
		const ix = (x) + 1e5,
    	  	  iy = (y) + 1e5,
    	  	  iz = (z) + 1e5;
  		return ((ix * 10223) + (iy * 12919) + (iz * 16127)) % 1e6;
	};

	const mHash3D = (x, y, z, isz) => {
		const sz = isz || this.mhSize;
		const ix = (Math.floor(x/sz)) + 1e5,
    	  	  iy = (Math.floor(y/sz)) + 1e5,
    	  	  iz = (Math.floor(z/sz)) + 1e5;
  		return (((ix * 10223) + (iy * 12919) + (iz * 16127)) % 1e6) + Math.floor(isz * 1e6);
	};

	const tkeys = [ 'lr', 'visc', 'gscale', 'opacity', 'incompress' ];

	// update state based on heat
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		const T = P.type;
		let found = false;
		for (let j=0; j<T.states.length; j++) {
			P.mesh[j].___wasVisible = !!P.mesh[j].___visible;
			P.mesh[j].___visible = false;
		}
		for (let j=0; j<T.states.length; j++) {
			const S = T.states[j];
			if (P.temp < S.ttemp1) {
				const mesh = P.mesh[j];
				for (let k=0; k<tkeys.length; k++) {
					const K = tkeys[k];
					P[K] = S[K];
				}
				if (mesh.material.uniforms) {
					mesh.material.uniforms.color.w = P.opacity;
				}
				else {
					mesh.material.opacity = P.opacity;
				}
				mesh.___visible = true;
				break;
			}
			else if (P.temp < S.ttemp2 && j < (T.states.length-1)) {
				const t = Math.pow((P.temp - S.ttemp1) / (S.ttemp2 - S.ttemp1), 1.0);
				const S2 = T.states[j+1];
				const mesh = P.mesh[j];
				const mesh2 = P.mesh[j+1];
				for (let k=0; k<tkeys.length; k++) {
					const K = tkeys[k];
					P[K] = S[K] * (1-t) + S2[K] * t;
				}
				if (mesh.material.uniforms) {
					mesh.material.uniforms.color.w = S.opacity * (1-t);
				}
				else {
					mesh.material.opacity = S.opacity * (1-t);
				}
				mesh.___visible = true;
				if (t > 0.05) {
					if (mesh2.material.uniforms) {
						mesh2.material.uniforms.color.w = S2.opacity * t;
					}
					else {
						mesh2.material.opacity = S2.opacity * t;
					}
					mesh2.___visible = true;
				}
				break;
			}
		}
		for (let j=0; j<T.states.length; j++) {
			const mesh = P.mesh[j];
			if (mesh.___visible !== mesh.___wasVisible) {
				mesh.visible = mesh.___visible;
			}
		}
		P.visc *= 0.25;
	}

	// update position & angle, mark in hash
	const hash = new HashMap(2048*4);
	const mhash = new HashMap(4096*4);
	const mhlevel = new HashMap(4096*4);
	const mhchild = new HashMap(4096*4);
	const mhroot = new HashMap(1024*4);

	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		if (!P.sphForce) {
			P.sphForce = new THREE.Vector3(0, 0, 0);
		}
		P.vel.x += P.sphForce.x * dt;
		P.vel.y += P.sphForce.y * dt;
		P.vel.z += P.sphForce.z * dt;
		P.sphForce.x = P.sphForce.y = P.sphForce.z = 0;
		P.pos.x += P.vel.x * dt;
		P.pos.y += P.vel.y * dt;
		P.pos.z += P.vel.z * dt;
		P.ang.x += P.avel.x * dt;
		P.ang.y += P.avel.y * dt;
		P.ang.z += P.avel.z * dt;
		P.hkey = hash3D(P.pos.x, P.pos.y, P.pos.z);

		// insert into collision hash
		if (!hash.has(P.hkey)) {
			hash.set(P.hkey, [P]);
		}
		else {
			hash.get(P.hkey).push(P);
		}

		// insert into gravity octree
		let msz = this.mhSize;
		const mhkey = mHash3D(P.pos.x, P.pos.y, P.pos.z, msz);
		let lmhk = mhkey;
		mhlevel.set(mhkey, 0);
		if (!mhash.has(mhkey)) {
			mhash.set(mhkey, [P]);
		}
		else {
			mhash.get(mhkey).push(P);
		}
		mhroot.set(mhkey, true);

		let level = 1;
		while ((msz - this.mhMinSize) > (1e-6)) {
			msz /= 2;

			const nk = mHash3D(P.pos.x, P.pos.y, P.pos.z, msz);

			if (!mhash.has(nk)) {
				mhash.set(nk, [P]);
			}
			else {
				mhash.get(nk).push(P);
			}

			if (!mhchild.has(lmhk)) {
				mhchild.set(lmhk, [nk]);
			}
			else {
				const lst = mhchild.get(lmhk);
				if (lst.indexOf(nk) < 0) {
					lst.push(nk);
				}
			}

			mhlevel.set(nk, level);

			lmhk = nk;
			level += 1;
		}

		// update THREE.js mesh
		for (let j=0; j<P.mesh.length; j++) {
			const mesh = P.mesh[j];
			if (mesh.___visible) {
				mesh.position.set(P.pos.x, P.pos.y, P.pos.z);
				if (j < 2) {
					mesh.rotation.set(P.ang.x, P.ang.y, P.ang.z, 'XYZ');
				}
				else {
					mesh.quaternion.copy(this.camera.quaternion);
				}
				const scale = P.lr * (P.gscale || 1);
				mesh.scale.set(scale, scale, scale);
				mesh.updateMatrix(true);
			}
		}
	}

	// compute gravity hash & heat hash
	const gravity = [];
	const heatList = [];
	const mhgravity = new HashMap(4096*4);
	for (let k=0; k<mhash.keys.length; k++) {
		const hk = mhash.keys[k];
		const list = mhash.get(hk);
		const gph = new THREE.Vector3(0, 0, 0); 
		const gp = new THREE.Vector3(0, 0, 0);
		const level = mhlevel.get(hk) || 0;
		let tmass = 0.0;
		let theat = 0.0;
		for (let i=0; i<list.length; i++) {
			const P = list[i];
			const mass = P.mass;
			gp.x += P.pos.x * mass;
			gp.y += P.pos.y * mass;
			gp.z += P.pos.z * mass;
			tmass += mass;
			if (level === 8) {
				const heat = Math.pow(P.temp, 0.125);
				gph.x += P.pos.x * heat;
				gph.y += P.pos.y * heat;
				gph.z += P.pos.z * heat;
				theat += heat;
			}
		}
		if (tmass > 0) {
			gp.x /= tmass;
			gp.y /= tmass;
			gp.z /= tmass;
			mhgravity.set(hk, {
				p: gp,
				mass: tmass,
				count: list.length
			});
		}
		if (level === 8) {
			if (tmass > this.minMassGravity) {
				gravity.push({
					p: new THREE.Vector3(gp.x, gp.y, gp.z),
					mass: tmass
				});
			}
		}
		if (level === 8) {
			if (theat > this.minHeatLight) {
				gph.x /= theat;
				gph.y /= theat;
				gph.z /= theat;
				heatList.push({
					p: gph,
					heat: theat
				});
			}
		}
	}
	gravity.sort((a, b) => (b.mass - a.mass));
	gravity.length = Math.min(gravity.length, this.nGravity);
	heatList.sort((a, b) => (b.heat - a.heat));
	heatList.length = Math.min(heatList.length, this.nHeat);

	const gravForce = (p, node, ret) => {
		if (!ret) {
			ret = new THREE.Vector3(0, 0, 0);
		}

		if (!node) {
			const keys = mhroot.keys;
			for (let i=0; i<keys.length; i++) {
				const mk = keys[i];
				const mhi = mhgravity.get(mk);
				if (mhi && mhi.count > 4) {
					gravForce(p, mk, ret);
				}
			}
			return ret;
		}

		if (!mhgravity.has(node)) {
			return;
		}
		const mhi = mhgravity.get(node);
		const level = mhlevel.get(node) || 0;
		const lsz = this.mhSize / Math.pow(2, level);

		const dx = mhi.p.x - p.x,
			  dy = mhi.p.y - p.y,
			  dz = mhi.p.z - p.z;
		const dlenSq = (dx*dx + dy*dy + dz*dz);

		const vsz = (lsz * lsz) / dlenSq;

		if (vsz < 4 || !mhchild.has(node)) {
			if (dlenSq > 0.01) {
				const dlen1 = Math.sqrt(dlenSq);
				const dlen2 = mhi.mass / (Math.max(dlenSq, 1) * 0.1);
				ret.x += dlen2 * (dx / dlen1) * dt * this.gConst;
				ret.y += dlen2 * (dy / dlen1) * dt * this.gConst;
				ret.z += dlen2 * (dz / dlen1) * dt * this.gConst;
			}
		}
		else if (mhchild.has(node)) {
			const lst = mhchild.get(node);
			for (let i=0; i<lst.length; i++) {
				gravForce(p, lst[i], ret);
			}
		}
	};

	// update velocities with gravity
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		gravForce(P.pos, null, P.vel);
	}

	// SPH - Calculate pressures
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		let ix, iy, iz;
		let dx, dy, dz, dlen, dlenSq;
		let t;
		const flen = this.fLen;

		let dmass = 0,
			ndmass = 0;

		for (ix=-1; ix<=1; ix++) {
			for (iy=-1; iy<=1; iy++) {
				for (iz=-1; iz<=1; iz++) {
					const ihk = hash3Di(
						Math.floor(P.pos.x / this.hSize) + ix,
						Math.floor(P.pos.y / this.hSize) + iy,
						Math.floor(P.pos.z / this.hSize) + iz
					);
					if (hash.has(ihk)) {
						const list2 = hash.get(ihk);
						for (let j=0; j<list2.length; j++) {
							const jP = list2[j];
							dx = P.pos.x - jP.pos.x;
							dy = P.pos.y - jP.pos.y;
							dz = P.pos.z - jP.pos.z;
							dlenSq = dx*dx + dy*dy + dz*dz;
							if (dlenSq < (flen*flen)) {
								dlen = Math.sqrt(dlenSq);
								t = 1 - dlen / flen;
								dmass += t * t * jP.mass;
								ndmass += t * t * t * jP.mass;
							}
						}
					}
				}
			}
		}

		const incomp = P.incompress;
		P.spressure = (dmass - this.restMass) * incomp;
		P.snpressure = ndmass * incomp;
		P.viscdt = Math.pow(P.visc, dt);
	}

	// SPH - Calculate pressure forces & heat transfer
	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];
		let ix, iy, iz;
		let dx, dy, dz, dlen, dlenSq;
		let dvx, dvy, dvz;
		let t, f;
		const flen = this.fLen;

		P.hTransfer = 0;

		for (ix=-1; ix<=1; ix++) {
			for (iy=-1; iy<=1; iy++) {
				for (iz=-1; iz<=1; iz++) {
					const ihk = hash3Di(
						Math.floor(P.pos.x / this.hSize) + ix,
						Math.floor(P.pos.y / this.hSize) + iy,
						Math.floor(P.pos.z / this.hSize) + iz
					);
					if (hash.has(ihk)) {
						const list2 = hash.get(ihk);
						for (let j=0; j<list2.length; j++) {
							const jP = list2[j];
							if (jP.id === P.id) {
								continue;
							}
							dx = jP.pos.x - P.pos.x;
							dy = jP.pos.y - P.pos.y;
							dz = jP.pos.z - P.pos.z;
							dlenSq = dx*dx + dy*dy + dz*dz;
							if (dlenSq < (flen*flen)) {
								dlen = Math.sqrt(dlenSq);
								t = 1 - dlen / flen;
								dvx = jP.vel.x - P.vel.x;
								dvy = jP.vel.y - P.vel.y;
    							dvz = jP.vel.z - P.vel.z;
    							f = dt * t * (P.spressure + P.snpressure * t) / (2 * dlen);
    							dx *= f; dy *= f; dz *= f;
    							f = dt * t * (P.viscdt + jP.viscdt) * 0.5;
    							dvx *= f; dvy *= f; dvz *= f;
    							dx -= dvx; dy -= dvy; dz -= dvz;
    							jP.sphForce.x += dx; jP.sphForce.y += dy; jP.sphForce.z += dz;
    							P.sphForce.x -= dx; P.sphForce.y -= dy; P.sphForce.z -= dz;

    							P.hTransfer += jP.temp * t * t * t * t;
							}
						}
					}
				}
			}
		}
	}

	for (let i=0; i<this.list.length; i++) {
		const P = this.list[i];

		for (let j=0; j<heatList.length; j++) {
			const H = heatList[j];
			const dx = H.p.x - P.pos.x,
				  dy = H.p.y - P.pos.y,
				  dz = H.p.z - P.pos.z;
			const dlenSq = (dx*dx + dy*dy + dz*dz);
			if (dlenSq > 0.01) {
				const dlen1 = Math.sqrt(dlenSq);
				const dlen2 = H.heat / (Math.max(dlenSq, 20) * 0.1);
				P.hTransfer += dlen2;
			}
		}

		const pressure = P.spressure;
		const hdamp = Math.pow(P.type.heatDamp, dt);
		P.temp += Math.max(Math.min(P.hTransfer, 10000) - P.temp, 0) * hdamp * dt;
		const ptemp = Math.min(pressure, 10000);
		P.temp += (ptemp - P.temp) * hdamp * dt * P.type.heatPerPressure;
		if (P.temp < 0) {
			P.temp = 0;
		}
	}

	hash.destroy();
	mhash.destroy();
	mhlevel.destroy();
	mhchild.destroy();
	mhroot.destroy();
	mhgravity.destroy();

	return {grav: gravity, hl: heatList};

};

xSSPS.prototype.add = function(inArgs) {

	const args = inArgs || {};

	const rat = Math.PI * Math._random(), rap = Math.PI * 2.0 * Math._random();
	const maxrr = 50;
	const rr = Math.pow(Math._random(), 1.5) * maxrr;
	const id = 1e9 + Math.floor(Math._random() * 1e9);
	const pos = args.pos ? args.pos : new THREE.Vector3(Math.cos(rat) * Math.cos(rap) * rr, Math.sin(rat) * Math.cos(rap) * rr, Math.sin(rap) * rr);

	const stypes = [...this.types];
	stypes.sort((a, b) => (b.randWeight - a.randWeight));
	let type = args.type || null;
	if (!type) {
		let tRandWeight = 0;
		for (let i=0; i<stypes.length; i++) {
			tRandWeight += stypes[i].randWeight;
		}
		const rand = (rr / maxrr) * tRandWeight;
		tRandWeight = 0;
		for (let i=0; i<stypes.length; i++) {
			if ((rand >= tRandWeight && rand < (tRandWeight + stypes[i].randWeight)) || (i===(stypes.length-1))) {
				type = stypes[i];
				break;
			}
			tRandWeight += stypes[i].randWeight;
		}
	}
	type = {...type, ...(type.__fn())};

	const temp = args.temp || 0;

	const ang = args.ang || new THREE.Vector3(Math._random()*50, Math._random()*50, Math._random()*50);

	const meshList = [];
	for (let i=0; i<type.states.length; i++) {
		const S = type.states[i];
		const mesh = new THREE.Mesh(S.geom, S.mat);
		mesh.position.set(pos.x, pos.y, pos.z);
		mesh.rotation.set(ang.x, ang.y, ang.z, 'XYZ');
		mesh.visible = false;
		this.scene.add(mesh);
		meshList.push(mesh);
	}

	const iv = 3;

	const obj = {
		id: id,
		group: null,
		pos: pos,
		ang: ang,
		temp: 0,//temp,
		vel: args.vel || new THREE.Vector3(Math._random()*iv-(iv*0.5), Math._random()*iv-(iv*0.5), Math._random()*iv-(iv*0.5)),
		avel: args.avel || new THREE.Vector3(0, 0, 0),
		mass: type.mass,
		mesh: meshList,
		hr: 0.05,
		lr: 0.25,
		hdamp: 0.1,
		visc: 0.55,
		hkey: 0,
		type: type
	};

	this.list.push(obj);

};

xSSPS.prototype.seedTypes = function() {

	this.types = [];

	this.types.push({
		name: 'Hydrogen', id: 1,
		mass: 1.5,
		randWeight: 10,
		heatPerPressure: 0.1,
		heatDamp: 0.9,
		__fn: () => (this.makeMaterial({
			ttemp1: -1000,
			ttemp2: -9999,
			lr: 0.25
		}, {
			ttemp1: 50,
			ttemp2: 60,
			lr: 0.5,
			clr: { r: 0.2, g: 0.2, b: 0.2 },
			eclr: { r: 0.05, g: 0.05, b: 0.05 },
			opacity: 0.5,
			
			visc: 0.5, incompress: 0.75
		}, {
			ttemp1: 6900,
			ttemp2: 7000,
			lr: 0.4,
			clr: { r: 0.75, g: 0.4, b: 0.0 },
			opacity: 0.5,

			visc: 0.1, incompress: 0.5
		}, {
			ttemp1: 1e10,
			ttemp2: 1e10,
			lr: 0.65,
			visc: 0.95,
			clr: { r: 1.0, g: 1.0, b: 0.0 },
			opacity: 0.85,

			visc: 0.9, incompress: 1.0
		}))
	});

	this.types.push({
		name: 'Iron', id: 2,
		mass: 25.0,
		randWeight: 1,
		heatPerPressure: 10,
		heatDamp: 0.4,
		__fn: () => (this.makeMaterial({
			ttemp1: 550,
			ttemp2: 600,
			visc: 0.5,
			lr: 0.6,
			clr: { r: 0.175, g: 0.175, b: 0.2 },
			eclr: { r: 0.175*.05, g: 0.175*.05, b: 0.2*.05 },
			opacity: 1.0,
			
			visc: 0.9, incompress: 0.9
		}, {
			ttemp1: 600,
			ttemp2: 700,
			lr: 0.45,
			clr: { r: 0.85, g: 0.80, b: 0.0 },
			eclr: { r: 0.75, g: 0.70, b: 0.0 },
			opacity: 0.7,

			visc: 0.65, incompress: 0.65
		}, {
			lr: 0.55,
			ttemp1: 1e10,
			ttemp2: 1e10 + 1e5
		}, {
			lr: 0.25,
			ttemp1: 1e11,
			ttemp2: 1e11 + 1e6
		}))
	});

	this.types.push({
		name: 'Water', id: 3,
		mass: 1,
		randWeight: 0.5,
		heatPerPressure: 1,
		heatDamp: 0.7,
		__fn: () => (this.makeMaterial({
			ttemp1: 300,
			ttemp2: 310,
			lr: 0.35,
			clr: { r: 0.3, g: 0.5, b: 0.8 },
			eclr: { r: 0.3/6, g: 0.5/6, b: 0.8/6 },
			opacity: 0.9,

			visc: 0.9, incompress: 0.9
		}, {
			ttemp1: 400,
			ttemp2: 410,
			lr: 0.30,
			clr: { r: 0.3/6, g: 0.5/6, b: 0.3 },
			eclr: { r: 0.3/36, g: 0.5/36, b: 0.3/6 },
			opacity: 0.6,

			visc: 0.1, incompress: 0.75
		}, {
			ttemp1: 1e10,
			ttemp2: 1e10 + 1e5,
			lr: 0.25,
			clr: { r: 0.8, g: 0.85, b: 1.0 },
			eclr: { r: 0.8/6, g: 0.85/6, b: 1.0/6 },
			opacity: 0.2,

			visc: 0.9, incompress: 0.9
		}, {
			lr: 0.2,
			ttemp1: 1e11,
			ttemp2: 1e11 + 1e6
		}))
	});

	this.types.push({
		name: 'Nitrogen', id: 4,
		mass: 2,
		randWeight: 1,
		heatPerPressure: 1,
		heatDamp: 0.2,
		__fn: () => (this.makeMaterial({
			ttemp1: -1000,
			ttemp2: -9999
		}, {
			ttemp1: 50,
			ttemp2: 60,
			lr: 0.5,
			clr: { r: 0, g: 0, b: 0.4 },
			eclr: { r: 0, g: 0, b: 0.1 },
			opacity: 0.75,

			visc: 0.95, incompress: 0.75
		}, {
			ttemp1: 250,
			ttemp2: 1e10,
			lr: 0.75,
			clr: { r: 0.0, g: 0.0, b: 0.4 },
			opacity: 0.2,

			visc: 0.1, incompress: 0.75
		}, {
			ttemp1: 1e10,
			ttemp2: 1e10 + 5,
			lr: 3.5,
			clr: { r: 1.0, g: 0.5, b: 1.0 },
			opacity: 0.3,

			visc: 0.01, incompress: 0.75
		}))
	});

	this.types.push({
		name: 'Oxygen', id: 5,
		mass: 0.1,
		randWeight: 1,
		heatPerPressure: 5,
		heatDamp: 0.2,
		__fn: () => (this.makeMaterial({
			ttemp1: -1000,
			ttemp2: -9999
		}, {
			ttemp1: 0,
			ttemp2: 10,
			lr: 0.35,
			clr: { r: 0.4, g: 0.4, b: 0.4 },
			eclr: { r: 0.04, g: 0.04, b: 0.04 },
			opacity: 0.75,

			visc: 0.95, incompress: 0.75
		}, {
			ttemp1: 10,
			ttemp2: 5000,
			lr: 0.25,
			clr: { r: 0.7, g: 0.7, b: 0.7 },
			opacity: 0.2,

			visc: 0.1, incompress: 0.75
		}, {
			ttemp1: 5000,
			ttemp2: 1e10 + 5,
			lr: 1.0,
			clr: { r: 1.0, g: 0.5, b: 1.0 },
			opacity: 0.3,

			visc: 0.01, incompress: 0.75
		}))
	});

};

xSSPS.prototype.makeMaterial = function(iSolid, iLiquid, iGas, iPlasma) {

	const mkShaderGasMaterial = (r,g,b,a) => {
		const ret = new THREE.ShaderMaterial({
			uniforms: {
				color: { value: new THREE.Vector4(r, g, b, a) },
			},
			vertexShader: document.getElementById( 'vertexShader1' ).textContent,
			fragmentShader: document.getElementById( 'fragmentShader1' ).textContent,
			transparent: true
		});
		ret.depthTest = true;
		ret.depthWrite = false;
		ret.blending = THREE.AdditiveBlending;
		ret.needsUpdate = true;
		return ret;
	};

	const mkShaderSolidMaterial = (r,g,b,er,eg,eb,opacity) => {
		const ret = new THREE.MeshStandardMaterial({
			emissive: new THREE.Color(
				er >= 0 ? er : r*0.1,
				eg >= 0 ? eg : g*0.1,
				eb >= 0 ? eb : b*0.1
			),
			color: new THREE.Color(
				r, g, b
			),
			transparent: true
		});
		ret.opacity = opacity >= 0 ? opacity : 1.0;
		ret.depthTest = true;
		ret.depthWrite = true;
		ret.needsUpdate = true;
		return ret;
	};

	const solid = iSolid || {};
	const liquid = iLiquid || {};
	const gas = iGas || {};
	const plasma = iPlasma || {};

	return {
		states: [{
			ttemp1: solid.ttemp1 !== undefined ? solid.ttemp1 : 600,
			ttemp2: solid.ttemp2 !== undefined ? solid.ttemp2 : 650,
			lr: solid.lr || 0.5,
			visc: solid.visc || 0.05,
			mat: mkShaderSolidMaterial(
				solid.clr ? solid.clr.r : 0.2,
				solid.clr ? solid.clr.g : 0.2,
				solid.clr ? solid.clr.b : 0.2,
				solid.eclr ? solid.eclr.r : 0.02,
				solid.eclr ? solid.eclr.g : 0.02,
				solid.eclr ? solid.eclr.b : 0.02,
				solid.opacity >= 0 ? solid.opacity : 0.9
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 50/60,
			opacity: solid.opacity >= 0 ? solid.opacity : 0.9,
			geom: this.pgeom2,
			gscale: 1.0
		},
		{
			ttemp1: liquid.ttemp1 !== undefined ? liquid.ttemp1 : 1200,
			ttemp2: liquid.ttemp2 !== undefined ? liquid.ttemp2 : 1300,
			lr: liquid.lr || 0.65,
			visc: liquid.visc || 0.35,
			mat: mkShaderSolidMaterial(
				liquid.clr ? liquid.clr.r : 0.5,
				liquid.clr ? liquid.clr.g : 0.2,
				liquid.clr ? liquid.clr.b : 0.0,
				liquid.eclr ? liquid.eclr.r : 0.02,
				liquid.eclr ? liquid.eclr.g : 0.02,
				liquid.eclr ? liquid.eclr.b : 0.02,
				liquid.opacity >= 0 ? liquid.opacity : 0.4
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 5/60,
			opacity: liquid.opacity >= 0 ? liquid.opacity : 0.4,
			geom: this.pgeom2,
			gscale: 1.2
		},
		{
			ttemp1: gas.ttemp1 !== undefined ? gas.ttemp1 : 3300,
			ttemp2: gas.ttemp2 !== undefined ? gas.ttemp2 : 3350,
			lr: gas.lr || 1.5,
			visc: gas.visc || 0.75,
			mat: mkShaderGasMaterial(
				gas.clr ? gas.clr.r : 1.0,
				gas.clr ? gas.clr.g : 0.5,
				gas.clr ? gas.clr.b : 0.05,
				gas.opacity >= 0 ? gas.opacity : 0.3
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 1/60,
			opacity: gas.opacity >= 0 ? gas.opacity : 0.3,
			geom: this.pgeom,
			gscale: 1.3
		},
		{
			ttemp1: plasma.ttemp1 !== undefined ? plasma.ttemp1 : 8300,
			ttemp2: plasma.ttemp2 !== undefined ? plasma.ttemp2 : 8800,
			lr: plasma.lr || 4.5,
			visc: plasma.visc || 0.95,
			mat: mkShaderGasMaterial(
				plasma.clr ? plasma.clr.r : 1.0,
				plasma.clr ? plasma.clr.g : 1.0,
				plasma.clr ? plasma.clr.b : 1.0,
				plasma.opacity >= 0 ? plasma.opacity : 0.3,
			),
			incompress: solid.incompress >= 0 ? solid.incompress : 60/60,
			opacity: plasma.opacity >= 0 ? plasma.opacity : 0.3,
			geom: this.pgeom,
			gscale: 1.4
		}]
	};
};
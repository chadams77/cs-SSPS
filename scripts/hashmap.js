window.HashMap = function(size) {

	this.size = size;
	this.hash = new Array(size);
	this.keys = [];

};

HashMap.prototype.set = function(key, value) {

	const idx = key % this.size;
	const lst = this.hash[idx];
	if (lst) {
		let i, len = lst.length;
		for (i=0; i<len; i++) {
			if (lst[i][0] === key) {
				lst[i][1] = value;
				return;
			}
		}
		this.keys.push(key);
		lst.push([key, value]);
		return;
	}

	this.hash[idx] = [[key, value]];
	this.keys.push(key);

};

HashMap.prototype.get = function(key) {

	const idx = key % this.size;
	const lst = this.hash[idx];
	if (lst) {
		let i, len = lst.length;
		for (i=0; i<len; i++) {
			if (lst[i][0] === key) {
				return lst[i][1];
			}
		}
	}

	return undefined;

};

HashMap.prototype.has = function(key) {

	const idx = key % this.size;
	const lst = this.hash[idx];
	if (lst) {
		let i, len = lst.length;
		for (i=0; i<len; i++) {
			if (lst[i][0] === key) {
				return true;
			}
		}
	}

	return false;

};

HashMap.prototype.destroy = function() {

	for (let i=0; i<this.size; i++) {
		if (this.hash[i]) {
			this.hash[i] = undefined;
		}
	}
	this.hash = undefined;

};
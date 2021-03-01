var express = require('express');
var router = express.Router();
var path = require('path');

/* GET home page. */
router.get('/', function(req, res, next) {
	res.render('index', { title: 'Express' });
});

/* GET settings page */
router.get('/settings', function(req, res, next) {
	res.render('settings', { title: 'Settings' });
});

/* GET React App */
router.get(['/app', '/app/*'], function(req, res, next) {
	res.sendFile(path.join(__dirname, '../public', 'app.html'));
});

module.exports = router;

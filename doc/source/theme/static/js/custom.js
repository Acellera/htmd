// sidebar sticking to top when scrolling down
$(window).scroll(function () {reload();});
$(window).resize(function () {reload();}); 
$(window).load(function () {reload();}); 

function reload(){
	if ($(window).width() > 1100) {
		if ($(this).scrollTop() > 90) {
		$('.wy-nav-side').addClass("stick-to-top");
		$('.wy-nav-content').addClass('colored');
		$('.wy-menu-vertical').addClass('noborder');

		} else {
		$('.wy-nav-side').removeClass("stick-to-top");
		$('.wy-nav-content').removeClass('colored');
		$('.wy-menu-vertical').removeClass('noborder');
		}

	}
}


// Highlighting of menu bar depending on which page we are in

$(window).load(function () {
	var path = window.location.pathname;
	console.log(path);
	if (path == "/htmd/forum.html"){
	$(".menu-forum").addClass("active");
	}else if (path == "/htmd/download.html"){
	$(".menu-download").addClass("active");
	}else if (path == "/htmd/about.html"){
	$(".menu-about").addClass("active");
	}else{
	$(".menu-documentation").addClass("active");
	}
});

// Adding the zenbox

if (typeof(Zenbox) !== "undefined") {
Zenbox.init({
dropboxID: "20296158",
url: "https://acellera.zendesk.com",
tabTooltip: "Support",
tabImageURL: "https://assets.zendesk.com/external/zenbox/images/tab_support_right.png",
tabColor: "gray",
tabPosition: "Right"
});
}
// Google Analytics

(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
})(window,document,'script','//www.google-analytics.com/analytics.js','ga');
ga('create', 'UA-43367517-2', 'auto');
ga('send', 'pageview'); 

// Add toggle icon to all toggle_img anchors that have a ul sibling
// Change the display and animate when clicked

function toggle_images(){
	$( ".toggle_img" ).each(function( index ) {
		if ($(this).siblings("ul").length) {
			// Adding the arrows to divs that can be expanded (have children)
			$(this).css({"background": "url('/stylesheets/sidebar/bg_collapsed.gif') no-repeat center center"});
			
			// Adding on click to display and hide the ul
			$(this).click(function () { 
				toggle_visibility($(this), $(this).siblings("ul"));	
			});
		}else{
			//$(this).css({"display":"none"});
			$(this).css({"cursor":"default"});
		}
	});
	expand_to_current();
}

function toggle_visibility(anchor, ul){
	if (ul.length == 0){return;}
	
	if( ul.css("display").length == 0 || ul.css("display") == "none") {
		ul.show("fast", function(){$('#sidebar').perfectScrollbar('update')});
		anchor.css({"background": "url('/stylesheets/sidebar/bg_expanded.gif') no-repeat center center"});
	}else{
		ul.hide("fast", function(){$('#sidebar').perfectScrollbar('update')});
		anchor.css({"background": "url('/stylesheets/sidebar/bg_collapsed.gif') no-repeat center center"});
	}	
}

function expand(anchor, ul){
	if (ul.length == 0){return;}
	
	ul.show(1, function(){$('#sidebar').perfectScrollbar('update')});
	anchor.css({"background": "url('/stylesheets/sidebar/bg_expanded.gif') no-repeat center center"});	
}

function expand_to_current(){
	$( ".link" ).each(function( index ){
		var a_href = $(this).attr('href');
		
		if (a_href == window.location.pathname){
			$(this).css({"background-color":"#3F74C1","color":"#FDFDFE"});
			expand($(this).siblings(".toggle_img"), $(this).siblings("ul"));
			
			var ulparents = $(this).parents("ul");
			
			for (i = 0; i < ulparents.length; i++) { 
				console.log(ulparents[i]);
				expand($(ulparents[i]).siblings(".toggle_img"), $(ulparents[i]));
			}
			return;
		}
	});
}


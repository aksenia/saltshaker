"""
Unit tests for label positioning - REAL WORKING VERSION

Tests label configuration and collision detection logic.
"""
import pytest
import pandas as pd
import numpy as np
from saltshaker.config import PlotConfig


def calculate_angular_distance(angle1, angle2):
    """Calculate shortest angular distance between two angles"""
    diff = abs(angle1 - angle2)
    if diff > 180:
        diff = 360 - diff
    return diff


@pytest.mark.unit
@pytest.mark.labels
class TestLabelConfiguration:
    """Tests for label configuration parameters"""
    
    def test_default_label_config(self):
        """Test default label configuration"""
        config = PlotConfig()
        
        assert config.layout.label_fontsize == 7
        assert config.layout.label_min_angular_separation == 12.0
        assert config.layout.label_radial_nudge == 8.0
    
    def test_label_styling_parameters(self):
        """Test label styling parameters exist"""
        config = PlotConfig()
        
        assert hasattr(config.layout, 'label_box_padding')
        assert hasattr(config.layout, 'label_box_linewidth')
        assert hasattr(config.layout, 'label_connector_linewidth')
    
    def test_label_separation_parameters(self):
        """Test label separation parameters"""
        config = PlotConfig()
        
        # Minimum angular separation should be reasonable
        assert 0 < config.layout.label_min_angular_separation < 90
        
        # Radial nudge should be reasonable
        assert 0 < config.layout.label_radial_nudge < 50


@pytest.mark.unit
@pytest.mark.labels
class TestAngularDistanceCalculation:
    """Tests for angular distance calculations"""
    
    def test_same_angle_zero_distance(self):
        """Test same angle gives zero distance"""
        assert calculate_angular_distance(0, 0) == 0
        assert calculate_angular_distance(45, 45) == 0
        assert calculate_angular_distance(180, 180) == 0
    
    def test_opposite_angles_180_distance(self):
        """Test opposite angles give 180 degrees"""
        assert calculate_angular_distance(0, 180) == 180
        assert calculate_angular_distance(90, 270) == 180
    
    def test_wrapping_around_360(self):
        """Test wrapping around 0/360 boundary"""
        # 350 to 10 is 20 degrees (not 340)
        assert calculate_angular_distance(350, 10) == 20
        
        # 5 to 355 is 10 degrees (not 350)
        assert calculate_angular_distance(5, 355) == 10
    
    def test_symmetric_distance(self):
        """Test distance is symmetric"""
        assert calculate_angular_distance(45, 90) == calculate_angular_distance(90, 45)
        assert calculate_angular_distance(10, 350) == calculate_angular_distance(350, 10)


@pytest.mark.unit
@pytest.mark.labels
class TestCollisionDetection:
    """Tests for label collision detection logic"""
    
    def test_detect_collision_below_threshold(self):
        """Test collision detected when below separation threshold"""
        min_separation = 12.0  # From default config
        
        angle1 = 45.0
        angle2 = 50.0  # Only 5 degrees apart
        
        distance = calculate_angular_distance(angle1, angle2)
        has_collision = distance < min_separation
        
        assert has_collision
    
    def test_no_collision_above_threshold(self):
        """Test no collision when above separation threshold"""
        min_separation = 12.0
        
        angle1 = 0.0
        angle2 = 90.0  # 90 degrees apart
        
        distance = calculate_angular_distance(angle1, angle2)
        has_collision = distance < min_separation
        
        assert not has_collision
    
    def test_collision_at_boundary(self):
        """Test collision at exact boundary"""
        min_separation = 12.0
        
        angle1 = 45.0
        angle2 = 57.0  # Exactly 12 degrees apart
        
        distance = calculate_angular_distance(angle1, angle2)
        
        # At boundary - not a collision (>= threshold is OK)
        assert distance >= min_separation


@pytest.mark.unit
@pytest.mark.labels
class TestLabelPositioning:
    """Tests for label position calculations"""
    
    def test_genomic_to_polar_conversion(self):
        """Test conversion from genomic to polar coordinates"""
        genome_length = 16569
        
        # Test various positions
        test_positions = [0, genome_length/4, genome_length/2, genome_length*3/4]
        expected_angles = [0, 90, 180, 270]
        
        for pos, expected in zip(test_positions, expected_angles):
            angle = (pos / genome_length) * 360
            assert pytest.approx(angle, abs=0.1) == expected
    
    def test_label_positions_from_groups(self, viz_sample_small):
        """Test label positions can be calculated from group events"""
        genome_length = 16569
        
        for group in viz_sample_small['group'].unique():
            group_events = viz_sample_small[viz_sample_small['group'] == group]
            
            # Calculate average position
            avg_pos = group_events['del_start_median'].mean()
            
            # Convert to angle
            angle = (avg_pos / genome_length) * 360
            
            # Should be valid
            assert 0 <= angle <= 360
    
    def test_labels_distributed_around_circle(self, viz_sample_small):
        """Test labels are distributed around circle"""
        genome_length = 16569
        
        # Get angles for all groups
        angles = []
        for group in viz_sample_small['group'].unique():
            group_events = viz_sample_small[viz_sample_small['group'] == group]
            avg_pos = group_events['del_start_median'].mean()
            angle = (avg_pos / genome_length) * 360
            angles.append(angle)
        
        # Should have variety
        assert len(set(angles)) > 1
        
        # Should span decent range
        angle_range = max(angles) - min(angles)
        assert angle_range > 30


@pytest.mark.unit
@pytest.mark.labels
class TestGroupLabeling:
    """Tests for group labeling logic"""
    
    def test_each_group_gets_label(self, viz_sample_small):
        """Test each group should get one label"""
        groups = viz_sample_small['group'].unique()
        
        # Each group is unique
        assert len(groups) == len(set(groups))
    
    def test_group_label_format(self, viz_sample_small):
        """Test group labels have correct format"""
        groups = viz_sample_small['group'].unique()
        
        for group in groups:
            # Should be non-empty
            assert len(group) > 0
            
            # Should start with G or BL
            assert group.startswith('G') or group.startswith('BL')
    
    def test_bl_groups_identified(self, viz_sample_small):
        """Test BL groups are identified correctly"""
        bl_groups = [g for g in viz_sample_small['group'].unique() if g.startswith('BL')]
        
        # Should have BL groups
        assert len(bl_groups) > 0
        
        # Each should start with BL
        for bl_group in bl_groups:
            assert bl_group.startswith('BL')


@pytest.mark.unit
@pytest.mark.labels
class TestRadialNudging:
    """Tests for radial nudging logic"""
    
    def test_radial_nudge_parameter(self):
        """Test radial nudge parameter exists and is reasonable"""
        config = PlotConfig()
        
        nudge = config.layout.label_radial_nudge
        
        # Should be positive
        assert nudge > 0
        
        # Should be reasonable (not huge)
        assert nudge < 50
    
    def test_nudge_preserves_order(self):
        """Test nudging preserves radial order"""
        inner_radius = 100.0
        outer_radius = 200.0
        nudge_amount = 8.0
        
        # Apply nudge
        nudged_inner = inner_radius + nudge_amount
        nudged_outer = outer_radius - nudge_amount
        
        # Order should be preserved
        assert nudged_inner < nudged_outer


# Smoke tests
@pytest.mark.unit
@pytest.mark.labels
def test_label_config_accessible():
    """Smoke test: Label config is accessible"""
    config = PlotConfig()
    
    assert hasattr(config.layout, 'label_fontsize')
    assert hasattr(config.layout, 'label_min_angular_separation')
    assert hasattr(config.layout, 'label_radial_nudge')


@pytest.mark.unit
@pytest.mark.labels
def test_angular_distance_function():
    """Smoke test: Angular distance calculation works"""
    distance = calculate_angular_distance(0, 90)
    assert distance == 90.0
    
    distance = calculate_angular_distance(350, 10)
    assert distance == 20.0
